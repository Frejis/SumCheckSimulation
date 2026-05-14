use ark_ff::Field;
use ark_poly::{Polynomial, SparseMultilinearExtension};
use ark_std::rand::{SeedableRng, rngs::StdRng};
use ark_std::test_rng;
use crate::gkr::layer::LayerReductionMessage;
use crate::structures::data_structures::SumCheckVerifier;
use crate::util::line_point;

pub struct StandardVerifier<F: Field> {
    random_points_chosen: Vec<F>,
    claimed_value: F,
    max_degree: usize,
    rng: StdRng,
}

impl<F: Field> StandardVerifier<F> {
    /// Function should only be used for testing as it relies on test_reng().
    pub fn handle_layer_reduction_message(&self, msg: LayerReductionMessage<F>, s_i_plus_1: usize) -> (Vec<F>, F) {
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

impl<F: Field> StandardVerifier<F> {
    pub fn new(max_degree: usize, claimed_value: F) -> Self {
        StandardVerifier {
            random_points_chosen: Vec::new(),
            claimed_value,
            max_degree,
            rng: StdRng::seed_from_u64(42),
        }
    }
}

impl<F: Field> SumCheckVerifier<F> for StandardVerifier<F> {
    fn verify_degree(&self, fx: &SparseMultilinearExtension<F>) -> bool {
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

    fn check_claimed_value(&self, gx: &SparseMultilinearExtension<F>) -> bool {
        let checked_claim: F = gx.evaluations.iter().map(|(_, &v)| v).sum();
        checked_claim == self.claimed_value
    }

    fn handle_round(&mut self, fx: &SparseMultilinearExtension<F>) -> F {
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

    fn final_check(&self) {
        todo!()
    }
}