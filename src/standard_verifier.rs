use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, Polynomial};
use crate::data_structures::{Prover, Verifier};

pub struct StandardVerifier<F: Field> {
    random_points_chosen: Vec<F>,
    claimed_value: F,
    max_degree: usize,
}

impl<F: Field> StandardVerifier<F> {
    pub fn new(max_degree: usize, claimed_value: F) -> Self {
        StandardVerifier {
            random_points_chosen: Vec::new(),
            claimed_value,
            max_degree,
        }
    }
}

impl<F: Field> Verifier<F> for StandardVerifier<F> {
    fn verify_degree(&self, fx: &DenseMultilinearExtension<F>) -> bool {
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

    fn check_claimed_value(&self, gx: &DenseMultilinearExtension<F>) -> bool {
        let checked_claim = gx.evaluations[0] + gx.evaluations[1];
        let res = checked_claim == self.claimed_value;
        println!("Claimed value {res} is correct.");
        res
    }

    fn handle_round(&mut self, fx: &DenseMultilinearExtension<F>) -> F {
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
}