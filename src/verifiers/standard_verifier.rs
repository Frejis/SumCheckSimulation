use ark_ff::Field;
use ark_poly::Polynomial;
use crate::gkr::gkr_round::GKRRound;
use crate::gkr::layer::LayerReductionMessage;
use crate::structures::data_structures::SumCheckVerifier;
use crate::util::{_line_point, interpolate_univariate};

pub struct StandardVerifier<F: Field> {
    random_points_chosen: Vec<F>,
    claimed_value: F,
    max_degree: usize,
    _gkr_round: GKRRound<F>,
    rng: ark_std::rand::rngs::StdRng,
}

impl<F: Field> StandardVerifier<F> {
    pub fn random_points_chosen(&self) -> &Vec<F> {
        &self.random_points_chosen
    }
}

impl<F: Field> StandardVerifier<F> {
    pub fn new(max_degree: usize, claimed_value: F, gkr_round: GKRRound<F>) -> Self {
        use ark_std::rand::SeedableRng;
        StandardVerifier {
            random_points_chosen: Vec::new(),
            claimed_value,
            max_degree,
            _gkr_round : gkr_round,
            rng: ark_std::rand::rngs::StdRng::from_entropy(),
        }
    }

    pub fn reduce_two_claims_to_one<R: rand::Rng>(
        &self,
        b_star: &[F],
        c_star: &[F],
        msg: &LayerReductionMessage<F>,
        rng: &mut R,
    ) -> (Vec<F>, F) {
        let q0 = msg.qt()[0];
        let q1 = msg.qt()[1];
        assert_eq!(q0, msg.z1(), "q(0) != z1");
        assert_eq!(q1, msg.z2(), "q(1) != z2");

        let r_star = F::rand(rng);
        let r_line_restriced = _line_point(b_star, c_star, r_star);
        let wlr = self._gkr_round.vi.evaluate(&r_line_restriced);
        let next_claim = interpolate_univariate(msg.qt(), r_star);

        // Check evaluation
        assert_eq!(next_claim, wlr, "q(r*) should equal W(line(r*))");

        // 1/|F| the above check is successful but eh, what can you do.
        (r_line_restriced, next_claim)
    }
}

impl<F: Field> SumCheckVerifier<F> for StandardVerifier<F> {
    fn verify_degree(&self, fx: &Vec<F>) -> bool {
        fx.len() <= self.max_degree + 1
    }

    /// Sample a random field element and store it in the transcript
    fn get_random_field_element(&mut self) -> F {
        // Use the verifier's persistent RNG so we get different values each time
        let rand_element = F::rand(&mut self.rng);
        self.random_points_chosen.push(rand_element);
        rand_element
    }

    fn check_claimed_value(&self, fx: &Vec<F>) -> bool {
        let q0 = fx[0];
        let q1 = fx[1];
        q0 + q1 == self.claimed_value
    }

    fn handle_round(&mut self, fx: &Vec<F>) -> F {
        if !self.verify_degree(fx) {
            panic!("Degree verification failed")
        }
        if !self.check_claimed_value(fx) {
            panic!("Claim check failed")
        }
        let r_i = self.get_random_field_element();
        // The new claim for the next round is g(r_i)
        self.claimed_value = interpolate_univariate(fx, r_i);
        r_i
    }

    fn set_claim(&mut self, claim: F) {
        self.claimed_value = claim;
    }

    fn final_check(&self) {
        todo!()
    }
}