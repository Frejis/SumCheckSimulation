use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, SparseMultilinearExtension};
use rand::Rng;
use crate::gkr::layer::Layer;
use crate::structures::circuit_structures::{GKRCircuit, Gate, GateType};

impl<F: Field> GKRCircuit<F> {
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
                // For now, only use Mul gates since the protocol doesn't handle mixed gates yet
                let typ = GateType::Mul;

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
        let mut mult_terms = Vec::<(usize, F)>::new();

        // Each gate corresponds to exactly one location in the hypercube.
        for (gate_idx, gate) in self.gates.iter().enumerate() {
            // Combined index over bits (x,b,c) for this wiring triple.
            // (x, b, c) -> (gate_idx, gate.left, gate.right)
            // LSB to MSB: x bits, then b bits, then c bits.
            let combined_index =
                gate_idx | (gate.left << k_x) | (gate.right << (k_x + k_child));
            match gate.typ {
                GateType::Add => add_terms.push((combined_index, F::one())),
                GateType::Mul => mult_terms.push((combined_index, F::one())),
            }
        }

        let total_vars = k_x + 2 * k_child;
        let add_ext = SparseMultilinearExtension::from_evaluations(total_vars, &add_terms);
        let mul_ext = SparseMultilinearExtension::from_evaluations(total_vars, &mult_terms);
        (add_ext, mul_ext)
    }

    /// The multilinear extension of gate values at this layer.
    pub fn value_extension(&self, k_x: usize) -> DenseMultilinearExtension<F> {
        DenseMultilinearExtension::from_evaluations_vec(k_x, self.values.clone())
    }
}