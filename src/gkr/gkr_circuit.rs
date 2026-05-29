use ark_ff::Field;
use ark_poly::SparseMultilinearExtension;
use ark_std::test_rng;
use rand::Rng;
use crate::gkr::gkr_driver::log2_pow2;
use crate::gkr::layer::{InputLayer, Layer};
use crate::gkr::predicates::{AddPredicate, MultPredicate};
use crate::structures::circuit_structures::{GKRCircuit, Gate, GateType};

impl<F: Field> GKRCircuit<F> {
    /// Generate a random layered circuit with the given sizes.
    /// Each element of `layer_sizes` must be a power of 2.
    /// Æast layer will be the input layer.
    pub fn random<R: Rng>(layer_sizes: &[usize], rng: &mut R) -> Self {
        assert!(layer_sizes.len() >= 2);
        let mut circuit_layer: Vec<Layer> = Vec::new();
        // This will be generated from bottom to top to easily generate the gates.
        for i in 0..layer_sizes.len() - 1 {
            let mut layer = Vec::new();
            for _ in 0..layer_sizes[i] {
                let gate_predicate: GateType = random_gate_type();
                let left_input: usize = random_number(layer_sizes[i +1], rng);
                let right_input: usize = random_number(layer_sizes[i +1], rng);
                let gate = Gate::new(left_input, right_input, gate_predicate);
                layer.push(gate);
            }
            circuit_layer.push(Layer::new(layer));
        }
        GKRCircuit::new(circuit_layer)
    }
}

fn random_number<R: Rng>(max: usize, rng: &mut R) -> usize {
    rng.gen_range(0..max)
}

fn random_gate_type() -> GateType {
    if test_rng().gen_bool(0.50f64) {
        GateType::Add
    } else {
        GateType::Mul
    }
}

/// This function initializes the predicates for the prover. This is a precomputation of all
/// predicates for the circuit.
pub fn compute_predicates<F: Field>(circuit: &GKRCircuit<F>, input: &InputLayer<F>) -> Vec<(AddPredicate<F>, MultPredicate<F>)> {
    let mut predicates = Vec::new();
    for i in 0..circuit.layers.len() {
        let mut add_terms = Vec::<(usize, F)>::new();
        let mut mul_terms = Vec::<(usize, F)>::new();

        let layer = &circuit.layers[i];
        // k_i denotes the address space, which is 2^(S_i).
        // Where S_i denotes the amount of gates in layer i.
        let s_i = log2_pow2(layer.gates.len());
        let next_s_i = log2_pow2(get_next_layer_address_space(input, i, circuit));

        // For the predicate we need to create an index of (g, b, c)
        // Where g is the gate index, b is the left child and c is the right child.
        // This is g is the first index then b and then c so it works when fixing variables.

        // TODO: Refactor below into a method that makes sense.
        for (gate_idx, gate) in layer.gates.iter().enumerate() {
            let left_index = gate.left << s_i;
            let right_index = gate.right << (s_i + next_s_i);
            let index: usize = gate_idx | left_index | right_index;
            match gate.predicate {
                GateType::Add => add_terms.push((index, F::one())),
                GateType::Mul => mul_terms.push((index, F::one())),
            }
        }

        let total_vars = s_i + 2 * next_s_i;
        let (add_pred, mult_pred) = create_predicate(&mut add_terms, &mut mul_terms, total_vars);
        predicates.push((add_pred,mult_pred))
    }
    predicates
}

fn create_predicate<F: Field>(add_terms: &mut Vec<(usize, F)>, mul_terms: &mut Vec<(usize, F)>, total_vars: usize) -> (AddPredicate<F>, MultPredicate<F>) {
    let add_sparse = SparseMultilinearExtension::from_evaluations(total_vars, &*add_terms);
    let add_pred: AddPredicate<F> = AddPredicate::new(add_sparse);
    let mult_sparse = SparseMultilinearExtension::from_evaluations(total_vars, &*mul_terms);
    let mult_pred: MultPredicate<F> = MultPredicate::new(mult_sparse);
    (add_pred, mult_pred)
}

fn get_next_layer_address_space<F: Field>(input: &InputLayer<F>, curr_layer: usize, circuit: &GKRCircuit<F>) -> usize {
    if curr_layer == circuit.layers.len() - 1 {
        input.values.len()
    } else {
        circuit.layers[curr_layer + 1].gates.len()
    }
}