use ark_ff::{Field, Zero};
use ark_poly::{DenseMultilinearExtension, Polynomial};
use ark_std::test_rng;
use crate::gkr::layer::InputLayer;
use crate::structures::circuit_structures::GKRCircuit;

/// This file implements a GKR verifier
/// The big difference from the sum-check verifier is that when evaluating the last round it will
/// Get two claimed values from the prover and restrict this to a line. Unless it is the last layer
/// Then it will check that the multilinear extension on the layer and check if the claim is true.

pub struct GKRVerifier<F: Field> {
    circuit: GKRCircuit<F>,
    input: InputLayer<F>,
    output_claim: DenseMultilinearExtension<F>,
}

impl<F: Field> GKRVerifier<F> {
    pub fn new(circuit: GKRCircuit<F>, input: InputLayer<F>) -> Self {
        Self { circuit, input, output_claim: DenseMultilinearExtension::zero() }
    }
}

impl<F: Field> GKRVerifier<F> {
    
    /// Panics if the evaluation of the input layer does not give the claimed evaluation
    /// When evaluated the mle of the input layer is evaluated in the "random" gate position. 
    /// This function panics == Verifier rejects.
    pub(crate) fn verify_final_claimed_value_point(&self, gate: Vec<F>, claimed_evaluation: F) {
        let input_mle = self.input.value_extension();
        assert_eq!(input_mle.evaluate(&gate), claimed_evaluation);
    }
}

impl<F: Field> GKRVerifier<F> {
    
    /// As i am only writing this for testing it is important to note that this function uses a
    /// test_rng, and should as such not be referenced for future production code.
    pub fn random_gate(&mut self, p0: &DenseMultilinearExtension<F>, gate_length: usize) -> Vec<F> {
        self.output_claim = p0.clone();
        let mut res = Vec::new();
        for _ in 0..gate_length {
            let rand_element = F::rand(&mut test_rng());
            res.push(rand_element);
        }
        res
    }
}
