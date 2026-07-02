//! The standard sum-check verifier.

use ark_ff::Field;
use ark_poly::univariate::SparsePolynomial;
use ark_poly::{
    DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension,
};
use ark_std::rand::rngs::StdRng;
use ark_std::rand::SeedableRng;
use ark_std::test_rng;

use crate::sumcheck::{LayerReductionMessage, SumCheckVerifier};
use crate::util::line_point;

/// Sum-check verifier that records the random points it hands out so the
/// final (or layer-reduction) check can be evaluated at them.
///
/// WARNING: only for simulation — the RNG is deterministically seeded and not
/// cryptographically secure.
pub struct StandardVerifier<F: Field> {
    random_points_chosen: Vec<F>,
    claimed_value: F,
    max_degree: usize,
    rng: StdRng,
}

impl<F: Field> StandardVerifier<F> {
    pub fn claimed_value(&self) -> &F {
        &self.claimed_value
    }

    pub fn random_points_chosen(&self) -> Vec<F> {
        self.random_points_chosen.clone()
    }

    /// Checks the prover's layer-reduction message and reduces the two claims
    /// `z1 = W(b*)`, `z2 = W(c*)` to a single claim at a random point on the
    /// line through `b*` and `c*`. Returns the next gate and its claimed value.
    ///
    /// Only for simulation: the point on the line is drawn from `test_rng`.
    pub fn handle_layer_reduction_message(
        &self,
        msg: &LayerReductionMessage<F>,
        s_i_plus_1: usize,
    ) -> (Vec<F>, F) {
        let q0 = msg.qt().evaluate(&F::zero());
        let q1 = msg.qt().evaluate(&F::one());
        assert_eq!(q0, msg.z1(), "q(0) != z1");
        assert_eq!(q1, msg.z2(), "q(1) != z2");

        let b_star = &self.random_points_chosen[0..s_i_plus_1];
        let c_star = &self.random_points_chosen[s_i_plus_1..2 * s_i_plus_1];

        // TODO: find a better way to draw this than a fresh test_rng.
        let r_star = F::rand(&mut test_rng());

        let next_gate = line_point(b_star, c_star, r_star);
        let next_claim = msg.qt().evaluate(&r_star);

        (next_gate, next_claim)
    }
}

impl<F: Field> SumCheckVerifier<F> for StandardVerifier<F> {
    fn new(max_degree: usize, claimed_value: F) -> Self {
        StandardVerifier {
            random_points_chosen: Vec::new(),
            claimed_value,
            max_degree,
            rng: StdRng::seed_from_u64(42),
        }
    }

    fn verify_degree(&self, fx: &SparsePolynomial<F>) -> bool {
        fx.degree() < self.max_degree
    }

    fn get_random_field_element(&mut self) -> F {
        let rand_element = F::rand(&mut self.rng);
        self.random_points_chosen.push(rand_element);
        rand_element
    }

    fn check_claimed_value(&self, gx: &SparsePolynomial<F>) -> bool {
        let checked_claim: F = gx.iter().rfold(F::zero(), |acc, (_, elem)| *elem + acc);
        checked_claim == self.claimed_value
    }

    fn handle_round(&mut self, fx: &SparsePolynomial<F>) -> F {
        if !self.verify_degree(fx) {
            panic!("round polynomial exceeds the degree bound");
        }
        if !self.check_claimed_value(fx) {
            panic!("round polynomial is inconsistent with the current claim");
        }
        self.get_random_field_element()
    }

    fn set_claim(&mut self, claim: F) {
        self.claimed_value = claim;
    }

    fn final_check(
        &self,
        gate: &[F],
        add_pred: &SparseMultilinearExtension<F>,
        mult_pred: &SparseMultilinearExtension<F>,
        mle: DenseMultilinearExtension<F>,
        vrf_func: Vec<(usize, F)>,
    ) {
        let add_fixed = add_pred.fix_variables(gate).evaluate(&self.random_points_chosen);
        let mult_fixed = mult_pred.fix_variables(gate).evaluate(&self.random_points_chosen);
        let half = self.random_points_chosen.len() / 2;
        let mle_left = mle.evaluate(&self.random_points_chosen[0..half].to_vec());
        let mle_right = mle.evaluate(&self.random_points_chosen[half..].to_vec());
        let add_term = add_fixed * (mle_left + mle_right);
        let mult_term = mult_fixed * (mle_left * mle_right);
        let checked_claim = SparseMultilinearExtension::from_evaluations(1, &vrf_func)
            .evaluate(&vec![*self.random_points_chosen.last().unwrap()]);
        assert_eq!(checked_claim, add_term + mult_term);
    }
}
