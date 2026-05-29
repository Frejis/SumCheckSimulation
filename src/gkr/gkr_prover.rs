use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, SparseMultilinearExtension};
use crate::gkr::gkr_driver::log2_pow2;
use crate::gkr::layer::InputLayer;
use crate::gkr::predicates::{AddPredicate, MultPredicate};
use crate::structures::circuit_structures::{EvaluatedGKRCircuit, GKRCircuit, GateType};

#[derive(Clone)]
pub struct GKRProver<F: Field> {
    circuit: GKRCircuit<F>,
    predicates: Vec<(AddPredicate<F>, MultPredicate<F>)>, // For each layer in the circuit.
    evaluated_circuit: EvaluatedGKRCircuit<F>,
    input: InputLayer<F>,
}

impl<F: Field> GKRProver<F> {
    pub fn set_predicates(&mut self, predicates: Vec<(AddPredicate<F>, MultPredicate<F>)>) {
        self.predicates = predicates;
    }
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

    pub fn input(&self) -> &InputLayer<F>{ &self.input }
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

    pub fn get_output_claim(&mut self) -> SparseMultilinearExtension<F> {
        let mut evaluations = Vec::new();
        for (idx, value) in self.evaluated_circuit.layers[0].values.iter().enumerate() {
            evaluations.push((idx, *value))
        }
        let vars = log2_pow2(self.evaluated_circuit.layers[0].values.len());
        SparseMultilinearExtension::from_evaluations(vars, &evaluations)
    }
}