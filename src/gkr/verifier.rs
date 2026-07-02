//! The GKR verifier.
//!
//! Unlike the standalone sum-check verifier, at the end of each layer's
//! sum-check the GKR verifier receives two claimed values from the prover and
//! restricts them to a line. Only at the input layer does it evaluate the
//! multilinear extension itself to check the final claim.

use ark_ff::Field;
use ark_poly::{MultilinearExtension, Polynomial};
use ark_std::test_rng;

use crate::gkr::{InputLayer, LayerPredicates};
use crate::sumcheck::{LayerReductionMessage, StandardVerifier};

pub struct GKRVerifier<F: Field> {
    input: InputLayer<F>,
    /// Per-layer wiring predicates, progressively fixed at the gates the
    /// verifier picks between layers.
    predicates: Vec<LayerPredicates<F>>,
}

impl<F: Field> GKRVerifier<F> {
    pub fn new(input: InputLayer<F>, predicates: Vec<LayerPredicates<F>>) -> Self {
        Self { input, predicates }
    }

    /// Picks a random gate label of the given length.
    ///
    /// Only for simulation: uses the deterministic `test_rng`, which is not
    /// cryptographically secure.
    pub fn random_gate(&self, gate_length: usize) -> Vec<F> {
        let mut rng = test_rng();
        (0..gate_length).map(|_| F::rand(&mut rng)).collect()
    }

    /// Fixes the predicates of `layer` at the gate the verifier picked, so
    /// they can be evaluated at the sum-check points later. Calls for the
    /// input layer (which has no predicates) are ignored.
    pub fn set_gate(&mut self, gate: &[F], layer: usize) {
        if layer == self.predicates.len() {
            return;
        }
        let add_fixed = self.predicates[layer].add.fix_variables(gate);
        let mult_fixed = self.predicates[layer].mult.fix_variables(gate);
        self.predicates[layer] = LayerPredicates { add: add_fixed, mult: mult_fixed };
    }

    /// Checks that the last sum-check claim of `layer` is consistent with the
    /// prover's claimed values `z1 = W(b*)` and `z2 = W(c*)`:
    /// `mult(b*,c*)·z1·z2 + add(b*,c*)·(z1 + z2)` must equal the running claim.
    /// Panicking == the verifier rejects.
    pub fn confirm_last_sum_check_msg(
        &self,
        msg: &LayerReductionMessage<F>,
        vrf: &StandardVerifier<F>,
        layer: usize,
    ) {
        let z1 = msg.z1();
        let z2 = msg.z2();
        let predicates = &self.predicates[layer];
        let b_c_star = vrf.random_points_chosen();
        let add_i = predicates.add.evaluate(&b_c_star);
        let mult_i = predicates.mult.evaluate(&b_c_star);
        let own_eval = z1 * z2 * mult_i + add_i * (z1 + z2);
        assert_eq!(own_eval, *vrf.claimed_value());
    }

    /// Evaluates the input layer's MLE at the final gate and checks it against
    /// the claimed evaluation. Panicking == the verifier rejects.
    pub fn verify_final_claimed_value_point(&self, gate: &[F], claimed_evaluation: F) {
        let input_mle = self.input.value_extension();
        assert_eq!(input_mle.evaluate(&gate.to_vec()), claimed_evaluation);
    }
}
