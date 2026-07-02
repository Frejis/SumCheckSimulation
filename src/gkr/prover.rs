//! The GKR prover.

use ark_ff::Field;
use ark_poly::SparseMultilinearExtension;

use crate::gkr::{EvaluatedGKRCircuit, GKRCircuit, InputLayer, LayerPredicates};
use crate::util::log2_pow2;

/// The GKR prover: evaluates the circuit once up front and keeps the
/// per-layer wiring predicates plus all gate values needed to answer the
/// verifier's sum-check queries.
pub struct GKRProver<F: Field> {
    predicates: Vec<LayerPredicates<F>>,
    evaluated_circuit: EvaluatedGKRCircuit<F>,
}

impl<F: Field> GKRProver<F> {
    pub fn new(
        circuit: &GKRCircuit<F>,
        input: &InputLayer<F>,
        predicates: Vec<LayerPredicates<F>>,
    ) -> Self {
        Self {
            predicates,
            evaluated_circuit: circuit.evaluate_circuit(input),
        }
    }

    pub fn predicates(&self) -> &[LayerPredicates<F>] {
        &self.predicates
    }

    pub fn evaluated_circuit(&self) -> &EvaluatedGKRCircuit<F> {
        &self.evaluated_circuit
    }

    /// The MLE of the output layer, sent to the verifier as the claimed output.
    pub fn output_claim(&self) -> SparseMultilinearExtension<F> {
        let output = &self.evaluated_circuit.layers[0].values;
        let evaluations: Vec<(usize, F)> =
            output.iter().enumerate().map(|(idx, value)| (idx, *value)).collect();
        SparseMultilinearExtension::from_evaluations(log2_pow2(output.len()), &evaluations)
    }
}
