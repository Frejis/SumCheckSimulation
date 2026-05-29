use ark_ff::{Field, Zero};
use ark_poly::{MultilinearExtension, Polynomial, SparseMultilinearExtension};
use ark_std::test_rng;
use crate::gkr::layer::{InputLayer, LayerReductionMessage};
use crate::gkr::predicates::{AddPredicate, MultPredicate};
use crate::structures::circuit_structures::GKRCircuit;
use crate::verifiers::standard_verifier::StandardVerifier;

/// This file implements a GKR verifier
/// The big difference from the sum-check verifier is that when evaluating the last round it will
/// Get two claimed values from the prover and restrict this to a line. Unless it is the last layer
/// Then it will check that the multilinear extension on the layer and check if the claim is true.

#[derive(Clone)]
pub struct GKRVerifier<F: Field> {
    circuit: GKRCircuit<F>,
    input: InputLayer<F>,
    output_claim: SparseMultilinearExtension<F>,
    pub(crate) predicate: Vec<(AddPredicate<F>, MultPredicate<F>)>,
    mi: F,
}

impl<F: Field> GKRVerifier<F> {
    pub fn set_predicate(&mut self, predicate: Vec<(AddPredicate<F>, MultPredicate<F>)>) {
        self.predicate = predicate;
    }
}

impl<F: Field> GKRVerifier<F> {
    pub fn new(circuit: GKRCircuit<F>, input: InputLayer<F>) -> Self {
        Self { circuit, input, output_claim: SparseMultilinearExtension::zero(), mi: F::zero(), predicate: Vec::new() }
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

    /// For a layer i it sets the predicate for the gate g picked by the verifier for both add and
    /// mult predicate.
    pub fn set_gate(&mut self, gate: &[F], layer: usize) {
        print!("Setting gate for layer {layer}.");
        if layer == self.predicate.len() {
            println!("Skipped input layer!");
            return;
        }
        print!("\n");
        let add_fixed = self.predicate[layer].0.pred.fix_variables(gate);
        let mult_fixed  = self.predicate[layer].1.pred.fix_variables(gate);
        self.predicate[layer].0 = AddPredicate::new(add_fixed);
        self.predicate[layer].1 = MultPredicate::new(mult_fixed);
    }
}

impl<F: Field> GKRVerifier<F> {
    
    /// As i am only writing this for testing it is important to note that this function uses a
    /// test_rng, and should as such not be referenced for future production code.
    /// Also updates the predicates for layer 0.
    pub fn random_gate(&mut self, out_claim: &SparseMultilinearExtension<F>, gate_length: usize) -> Vec<F> {
        self.output_claim = out_claim.clone();
        let mut res = Vec::new();
        for _ in 0..gate_length {
            let rand_element = F::rand(&mut test_rng());
            res.push(rand_element);
        }
        res
    }

    pub fn check_next_layer_claim_is_sum(&self, sum: F) {
        assert_eq!(sum, self.mi)
    }

    pub fn set_next_layer_claim(&mut self, mi: F) {
        self.mi = mi
    }

    pub fn confirm_last_sum_check_msg(&self, msg: &LayerReductionMessage<F>, vrf: &StandardVerifier<F>, layer: usize) {
        let z1 = msg.z1();
        let z2 = msg.z2();
        let (add_i, mult_i) =  &self.predicate[layer];
        let b_c_star = vrf.random_points_chosen();
        let add_i = add_i.pred.evaluate(&b_c_star);

        let mult_i = mult_i.pred.evaluate(&b_c_star);
        let own_eval = z1*z2*mult_i + add_i*(z1 + z2);
        assert_eq!(own_eval, *vrf.claimed_value())
    }
}
