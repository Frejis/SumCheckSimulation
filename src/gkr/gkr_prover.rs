use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, SparseMultilinearExtension};
use crate::gkr::gkr_driver::log2_pow2;
use crate::gkr::layer::InputLayer;
use crate::gkr::predicates::{AddPredicate, MultPredicate};
use crate::structures::circuit_structures::{EvaluatedGKRCircuit, GKRCircuit, GateType};

pub struct GKRProver<F: Field> {
    circuit: GKRCircuit<F>,
    predicates: Vec<(AddPredicate<F>, MultPredicate<F>)>, // For each layer in the circuit.
    evaluated_circuit: EvaluatedGKRCircuit<F>,
    input: InputLayer<F>,
}

impl<F: Field> GKRProver<F> {
    pub fn evaluated_circuit(circuit: &GKRCircuit<F>, input: &InputLayer<F>) -> EvaluatedGKRCircuit<F> {
        circuit.evaluate_circuit(input)
    }

    pub fn circuit(&self) -> &GKRCircuit<F> {
        &self.circuit
    }

    pub fn eval_circuit(&self) -> &EvaluatedGKRCircuit<F> {&self.evaluated_circuit}

    pub fn predicates(&self) -> &Vec<(AddPredicate<F>, MultPredicate<F>)> {
        &self.predicates
    }
}

impl<F: Field> GKRProver<F> {
    pub fn new(circuit: GKRCircuit<F>, input: InputLayer<F>) -> Self {
        let evaluated_circuit = GKRProver::evaluated_circuit(&circuit, &input);
        Self {
            circuit,
            predicates: Vec::new(),
            evaluated_circuit,
            input,
        }
    }

    /// This function initializes the predicates for the prover. This is a precomputation of all
    /// predicates for the circuit.
    pub fn compute_predicates(&mut self) {
        for i in 0..self.circuit.layers.len() {
            let mut add_terms = Vec::<(usize, F)>::new();
            let mut mul_terms = Vec::<(usize, F)>::new();

            let layer = &self.circuit.layers[i];
            // k_i denotes the address space, which is 2^(S_i).
            // Where S_i denotes the amount of gates in layer i.
            let s_i = log2_pow2(layer.gates.len());
            let next_s_i = log2_pow2(self.get_next_layer_address_space(i));

            // For the predicate we need to create an index of (g, b, c)
            // Where g is the gate index, b is the left child and c is the right child.
            // This is g is the first index then b and then c so it works when fixing variables.

            // TODO: Refactor below into a method that makes sense.
            for (gate_idx, gate) in layer.gates.iter().enumerate() {
                let left_index = gate.left << s_i;
                let right_index = gate.right << s_i + next_s_i;
                let index: usize = gate_idx | left_index | right_index;
                match gate.predicate {
                    GateType::Add => add_terms.push((index, F::one())),
                    GateType::Mul => mul_terms.push((index, F::one())),
                }
            }

            let total_vars = s_i + 2 * next_s_i;
            let (add_pred, mult_pred) = Self::create_predicate(&mut add_terms, &mut mul_terms, total_vars);
            self.predicates.push((add_pred,mult_pred))
        }
    }

    fn create_predicate(add_terms: &mut Vec<(usize, F)>, mul_terms: &mut Vec<(usize, F)>, total_vars: usize) -> (AddPredicate<F>, MultPredicate<F>) {
        let add_sparse = SparseMultilinearExtension::from_evaluations(total_vars, &*add_terms);
        let add_pred: AddPredicate<F> = AddPredicate::new(add_sparse);
        let mult_sparse = SparseMultilinearExtension::from_evaluations(total_vars, &*mul_terms);
        let mult_pred: MultPredicate<F> = MultPredicate::new(mult_sparse);
        (add_pred, mult_pred)
    }

    fn get_next_layer_address_space(&self, i: usize) -> usize {
        if i == self.circuit.layers.len() - 1 {
            self.input.values.len()
        } else {
            self.circuit.layers[i + 1].gates.len()
        }
    }

    pub fn get_output_claim(&mut self) -> SparseMultilinearExtension<F> {
        let mut evaluations = Vec::new();
        for (idx, value) in self.evaluated_circuit.layers[0].values.iter().enumerate() {
            evaluations.push((idx, *value))
        }
        let vars = log2_pow2(self.evaluated_circuit.layers[0].values.len());
        SparseMultilinearExtension::from_evaluations(vars, &evaluations)
    }
}