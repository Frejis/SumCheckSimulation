use ark_ff::Field;
use ark_poly::{Polynomial, SparseMultilinearExtension};
use ark_std::test_rng;
use crate::gkr::gkr_round::GKRRound;
use crate::gkr::layer::LayerReductionMessage;
use crate::structures::data_structures::SumCheckVerifier;
use crate::util::_line_point;

pub struct StandardVerifier<F: Field> {
    random_points_chosen: Vec<F>,
    claimed_value: F,
    max_degree: usize,
    _gkr_round: GKRRound<F>,
}

impl<F: Field> StandardVerifier<F> {
    pub fn new(max_degree: usize, claimed_value: F, gkr_round: GKRRound<F>) -> Self {
        StandardVerifier {
            random_points_chosen: Vec::new(),
            claimed_value,
            max_degree,
            _gkr_round : gkr_round,
        }
    }

    pub fn reduce_two_claims_to_one<R: rand::Rng>(
        &self,
        b_star: &[F],
        c_star: &[F],
        msg: &LayerReductionMessage<F>,
        rng: &mut R,
    ) -> (Vec<F>, F) {
        let q0 = msg.qt().evaluate(&[F::zero()].to_vec());
        let q1 = msg.qt().evaluate(&[F::one()].to_vec());
        assert_eq!(q0, msg.z1(), "q(0) != z1");
        assert_eq!(q1, msg.z2(), "q(1) != z2");

        let r_star = F::rand(&mut test_rng()); // TODO find a better way to do this.
        let r_line_restriced = _line_point(b_star, c_star, r_star);
        let wlr = self._gkr_round.vi.evaluate(&r_line_restriced);

        let next_claim = msg.qt().evaluate(&[r_star].to_vec());

        // Check evaluation
        assert_eq!(next_claim, wlr);

        // 1/|F| the above check is successful but eh, what can you do.
        (r_line_restriced, next_claim)
    }
}

impl<F: Field> SumCheckVerifier<F> for StandardVerifier<F> {
    fn verify_degree(&self, fx: &SparseMultilinearExtension<F>) -> bool {
        fx.degree() < self.max_degree
    }

    /// WARNING THIS SHOULD ONLY BE USED FOR TESTING SINCE TEST_RNG IS NOT CRYPTOGRAPHICALLY SECURE
    /// This is also just for simulating, so for now it will be fine to give an idea.
    fn get_random_field_element(&mut self) -> F {
        let mut rng = ark_std::test_rng();
        // Let's sample uniformly random field elements:
        let rand_element = F::rand(&mut rng);
        self.random_points_chosen.push(rand_element);
        rand_element
    }

    fn check_claimed_value(&self, gx: &SparseMultilinearExtension<F>) -> bool {
        let checked_claim: F = gx.evaluations.iter().map(|(_, &v)| v).sum();
        let res = checked_claim == self.claimed_value;
        res
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