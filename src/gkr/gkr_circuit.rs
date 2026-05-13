use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, SparseMultilinearExtension};
use ark_std::test_rng;
use rand::Rng;
use crate::gkr::layer::Layer;
use crate::structures::circuit_structures::{GKRCircuit, Gate, GateType};

impl<F: Field> GKRCircuit<F> {
    /// Generate a random layered circuit with the given sizes.
    /// Each element of `layer_sizes` must be a power of 2.
    /// Æast layer will be the input layer.
    /// Should only be used for testing as it uses test_rng()
    pub fn random<R: Rng>(layer_sizes: &[usize], rng: &mut R) -> Self {
        assert!(layer_sizes.len() >= 2);
        let mut circuit_layer: Vec<Layer> = Vec::new();
        // This will be generated from bottom to top to easily generate the gates.
        for i in 0..layer_sizes.len() - 1 {
            let mut layer = Vec::new();
            for j in 0..layer_sizes[i] {
                let gate_predicate: GateType = random_gate_type();
                let left_input: usize = random_number(layer_sizes[i +1]);
                let right_input: usize = random_number(layer_sizes[i +1]);
                let gate = Gate::new(left_input, right_input, gate_predicate);
                layer.push(gate);
            }
            circuit_layer.push(Layer::new(layer));
        }
        GKRCircuit::new(circuit_layer)
    }
}

fn random_number(max: usize) -> usize {
    test_rng().gen_range(0..max)
}

fn random_gate_type() -> GateType {
    if test_rng().gen_bool(0.50f64) {
        GateType::Add
    } else {
        GateType::Mul
    }
}