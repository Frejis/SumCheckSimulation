use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, SparseMultilinearExtension};
use rand::Rng;

/// Gate type: add or multiply child outputs.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum GateType {
    Add,
    Mul,
}

/// A gate at one layer, referencing children at the next layer.
#[derive(Clone)]
pub struct Gate {
    left: usize,
    right: usize,
    typ: GateType,
}

/// A layered arithmetic circuit.
pub struct GkrCircuit<F: Field> {
    pub layers: Vec<Layer<F>>,
}

pub struct Layer<F: Field> {
    /// Gate wiring for this layer (size = 2^{k_i}).
    pub gates: Vec<Gate>,
    /// Gate values at this layer.
    pub values: Vec<F>,
}

impl<F: Field> GkrCircuit<F> {
    /// Generate a random layered circuit with the given sizes.
    /// Each element of `layer_sizes` must be a power of 2.
    pub fn random<R: Rng>(layer_sizes: &[usize], rng: &mut R) -> Self {
        assert!(layer_sizes.len() >= 2);

        // Build values for input (last) layer randomly.
        let last_idx = layer_sizes.len() - 1;
        let mut layers = Vec::with_capacity(layer_sizes.len());
        layers.push(Layer {
            gates: Vec::new(),
            values: (0..layer_sizes[last_idx]).map(|_| F::rand(rng)).collect(),
        });

        // Build remaining layers bottom-up.
        for i in (0..last_idx).rev() {
            let next_size = layer_sizes[i + 1];
            let mut gates = Vec::with_capacity(layer_sizes[i]);
            let mut values = Vec::with_capacity(layer_sizes[i]);

            for _ in 0..layer_sizes[i] {
                let left = rng.gen_range(0..next_size);
                let right = rng.gen_range(0..next_size);
                let typ = if rng.r#gen::<bool>() {
                    GateType::Add
                } else {
                    GateType::Mul
                };

                // Compute this gate’s value from next layer’s values.
                let val = match typ {
                    GateType::Add => layers.last().unwrap().values[left]
                        + layers.last().unwrap().values[right],
                    GateType::Mul => layers.last().unwrap().values[left]
                        * layers.last().unwrap().values[right],
                };
                gates.push(Gate { left, right, typ });
                values.push(val);
            }

            layers.push(Layer { gates, values });
        }

        // Reverse so layer 0 is the output layer.
        layers.reverse();
        Self { layers }
    }
}

impl<F: Field> Layer<F> {
    /// Build the wiring predicate extensions for this layer.
    /// |x| = k_x = log2(#gates), |b|=|c|=k_child = log2(#child gates).
    pub fn wiring_predicates(
        &self,
        k_x: usize,
        k_child: usize,
    ) -> (SparseMultilinearExtension<F>, SparseMultilinearExtension<F>) {
        let mut add_terms = Vec::<(usize, F)>::new();
        let mut mul_terms = Vec::<(usize, F)>::new();

        // Each gate corresponds to exactly one location in the hypercube.
        for (gate_idx, gate) in self.gates.iter().enumerate() {
            // Combined index over bits (x,b,c) for this wiring triple
            let combined_index =
                (gate_idx << (2 * k_child)) | (gate.left << k_child) | gate.right;
            match gate.typ {
                GateType::Add => add_terms.push((combined_index, F::one())),
                GateType::Mul => mul_terms.push((combined_index, F::one())),
            }
        }

        let total_vars = k_x + 2 * k_child;
        let add_ext = SparseMultilinearExtension::from_evaluations(total_vars, &add_terms);
        let mul_ext = SparseMultilinearExtension::from_evaluations(total_vars, &mul_terms);
        (add_ext, mul_ext)
    }

    /// The multilinear extension of gate values at this layer.
    pub fn value_extension(&self, k_x: usize) -> DenseMultilinearExtension<F> {
        DenseMultilinearExtension::from_evaluations_vec(k_x, self.values.clone())
    }
}