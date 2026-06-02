use std::ops::Mul;
use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use ark_poly::univariate::SparsePolynomial;
use ark_std::rand::{SeedableRng, rngs::StdRng};
use ark_std::test_rng;
use crate::gkr::layer::LayerReductionMessage;
use crate::gkr::predicates::{AddPredicate, MultPredicate};
use crate::structures::data_structures::SumCheckVerifier;
use crate::util::line_point;

pub struct StandardVerifier<F: Field> {
    random_points_chosen: Vec<F>,
    claimed_value: F,
    max_degree: usize,
    rng: StdRng,
}

impl<F: Field> StandardVerifier<F> {
    pub(crate) fn claimed_value(&self) -> &F {
        &self.claimed_value
    }
}

impl<F: Field> StandardVerifier<F> {
    /// Function should only be used for testing as it relies on test_reng().
    pub fn handle_layer_reduction_message(&self, msg: &LayerReductionMessage<F>, s_i_plus_1: usize) -> (Vec<F>, F) {
        let q0 = msg.qt().evaluate(&F::zero());
        let q1 = msg.qt().evaluate(&F::one());
        assert_eq!(q0, msg.z1(), "q(0) != z1");
        assert_eq!(q1, msg.z2(), "q(1) != z2");

        let b_star = self.random_points_chosen[0..s_i_plus_1].to_vec();
        let c_star = self.random_points_chosen[s_i_plus_1..2*s_i_plus_1].to_vec();

        let r_star = F::rand(&mut test_rng()); // TODO find a better way to do this.

        let next_gate = line_point(&b_star, &c_star, r_star);

        let next_claim = msg.qt().evaluate(&r_star);

        (next_gate, next_claim)
    }
}

impl<F: Field> StandardVerifier<F> {
    pub(crate) fn random_points_chosen(&self) -> Vec<F> {
        self.random_points_chosen.clone()
    }
}

impl<F: Field> SumCheckVerifier<F> for StandardVerifier<F> {
    fn verify_degree(&self, fx: &SparsePolynomial<F>) -> bool {
        fx.degree() < self.max_degree
    }

    /// WARNING THIS SHOULD ONLY BE USED FOR TESTING SINCE TEST_RNG IS NOT CRYPTOGRAPHICALLY SECURE
    /// This is also just for simulating, so for now it will be fine to give an idea.
    fn get_random_field_element(&mut self) -> F {
        // Let's sample uniformly random field elements:
        let rand_element = F::rand(&mut self.rng);
        self.random_points_chosen.push(rand_element);
        rand_element
    }

    fn check_claimed_value(&self, gx: &SparsePolynomial<F>) -> bool {
        let checked_claim: F = gx.iter().rfold(F::zero(), |acc, (_, elem)| {*elem + acc});
        checked_claim == self.claimed_value
    }

    fn handle_round(&mut self, fx: &SparsePolynomial<F>) -> F {
        if !self.verify_degree(fx) {
            panic!()
        }
        if !self.check_claimed_value(fx) {
            panic!()
        }
        self.get_random_field_element()
    }

    fn set_claim(&mut self, claim: F) {
        self.claimed_value = claim;
    }

    fn final_check(&self, gate: &[F], add_pred: &SparseMultilinearExtension<F>, mult_pred: &SparseMultilinearExtension<F>, mle: DenseMultilinearExtension<F>, vrf_func: SparsePolynomial<F>) {
        println!("numvars: {:?}, input: {:?}", mle.num_vars, self.random_points_chosen.len());
        let add_fixed = add_pred.fix_variables(gate).evaluate(&self.random_points_chosen);
        let mult_fixed = mult_pred.fix_variables(gate).evaluate(&self.random_points_chosen);
        let mle_left = mle.evaluate(&self.random_points_chosen[0..self.random_points_chosen.len()/2].to_vec());
        let vj = &self.random_points_chosen[self.random_points_chosen.len()/2..self.random_points_chosen.len()];
        let mle_right = mle.evaluate(&vj.to_vec());
        let add_term = add_fixed * mle_left + add_fixed + mle_right;
        let mult_term = mult_fixed * (mle_left * mle_right);
        let checked_claim: F = vrf_func.iter().rfold(F::zero(), |acc, (_, elem)| {*elem + acc});
        assert_eq!(checked_claim, add_term + mult_term)
    }

    fn new(max_degree: usize, claimed_value: F) -> Self {
        StandardVerifier {
            random_points_chosen: Vec::new(),
            claimed_value,
            max_degree,
            rng: StdRng::seed_from_u64(42),
        }
    }
}