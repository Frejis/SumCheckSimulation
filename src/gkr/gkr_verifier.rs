use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, SparseMultilinearExtension};
use ark_std::test_rng;
use crate::gkr::layer::InputLayer;
use crate::gkr::predicates::{AddPredicate, MultPredicate};
use crate::structures::circuit_structures::{EvaluatedGKRCircuit, GKRCircuit};
use crate::structures::data_structures::{SumCheckProver, SumCheckVerifier};

/// This file implements a GKR verifier
/// The big difference from the sum-check verifier is that when evaluating the last round it will
/// Get two claimed values from the prover and restrict this to a line. Unless it is the last layer
/// Then it will check that the multilinear extension on the layer and check if the claim is true.

pub struct GKRVerifier<F: Field> {
    circuit: GKRCircuit<F>,
    input: InputLayer<F>,
    claimed_circuit_eval: DenseMultilinearExtension<F>,
    round: u32,
    output_claim: DenseMultilinearExtension<F>,
}

impl<F: Field> GKRVerifier<F> {
    
    /// As i am only writing this for testing it is important to note that this function uses a
    /// test_rng, and should as such not be referenced for future production code.
    pub fn random_gate(&mut self, p0: &DenseMultilinearExtension<F>, gate_length: usize) -> Vec<F> {
        self.output_claim = p0.clone();
        let mut res = Vec::new();
        for i in 0..gate_length {
            let rand_element = F::rand(&mut test_rng());
            res.push(rand_element);
        }
        res
    }
}
