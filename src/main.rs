use ark_bls12_381::Fr;
use ark_std::test_rng;
use crate::data_structures::{Prover, Verifier};
use crate::naive_sum_check::NaiveProver;
use crate::standard_verifier::StandardVerifier;
use crate::util::{random_gate, random_gkr_instance};

mod data_structures;
mod util;
pub mod naive_sum_check;
mod standard_verifier;
mod fast_prover;

fn main() {
    println!("Hello, world!");

    // Just testing things out. First create a random GKR instance with Fr
    let variables = 10;
    let mut rng = test_rng();
    let (mult, vi, vj) = random_gkr_instance::<Fr, _>(10, &mut rng);
    let rand_gate = random_gate(variables);

    let prover = NaiveProver::new(mult, vi, vj, rand_gate);

    // I have no clue about the max degree atm so i just say fuck it we ball hahaha
    let mut verifier = StandardVerifier::new(variables*400, prover.compute_sum());
    let verifier_func = prover.get_verifier_function();
    if verifier.check_claimed_value(&verifier_func) {
        println!("Claimed value is correct.");
    }
}

mod generic_tests {
    use ark_bls12_381::Fr;
    use ark_std::{test_rng};
    use crate::data_structures::{Prover, Verifier};
    use crate::naive_sum_check::NaiveProver;
    use crate::standard_verifier::StandardVerifier;
    use crate::util;
    use crate::util::{random_gate};

    #[test]
    fn test_verifier_first_round() {
        // 7 variables 3 seconds.... 8 variables 23!!! seconds!??! :OOO whaaa

        let variables = 7;
        let fixed_gate = random_gate(variables);

        let mut rng = test_rng();
        let (mult, vi, vj) = util::random_gkr_instance(variables, &mut rng);
        let prover: NaiveProver<Fr> = NaiveProver::new(mult, vi, vj, Vec::from(fixed_gate));
        // Now we test the g_func gives what we expect
        let verifier = StandardVerifier::new(variables, prover.compute_sum());

        let verifier_func = prover.get_verifier_function();


        assert!(verifier.check_claimed_value(&verifier_func));
    }

    #[test]
    fn simulate_two_rounds_naive() {
        let variables = 7;
        let fixed_gate = random_gate(variables);

        let mut rng = test_rng();
        let (mult, vi, vj) = util::random_gkr_instance(variables, &mut rng);
        let mut prover: NaiveProver<Fr> = NaiveProver::new(mult, vi, vj, Vec::from(fixed_gate));
        // Now we test the g_func gives what we expect
        let mut verifier = StandardVerifier::new(variables, prover.compute_sum());

        let verifier_func = prover.get_verifier_function();

        assert!(verifier.check_claimed_value(&verifier_func));
        let rand_var = verifier.get_random_field_element();

        prover.fix_variable(rand_var);

        let new_claim = prover.compute_sum();
        verifier.set_claim(new_claim);
        let snd_verifier_func = prover.get_verifier_function();
        assert!(verifier.check_claimed_value(&snd_verifier_func));
    }
}