use ark_bls12_381::Fr;
use crate::data_structures::{GKRRound, Prover, Verifier};
use crate::naive_sum_check::NaiveProver;
use crate::standard_verifier::StandardVerifier;
use crate::util::random_gate;

mod data_structures;
mod util;
pub mod naive_sum_check;
mod standard_verifier;
mod fast_prover;
pub mod circuit_structures;

fn main() {
    println!("Hello, world!");

    // Just testing things out. First create a random GKR instance with Fr
    let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
    let random_gate = random_gate(gkr_round.gate_labes());
    let mut prover = NaiveProver::new(gkr_round.clone(), &random_gate);
    // Now we test the g_func gives what we expect
    let verifier = StandardVerifier::new(3, prover.compute_sum(), gkr_round);
    let verifier_func = prover.get_verifier_function();
    if verifier.check_claimed_value(&verifier_func) {
        println!("Claimed value is correct.");
    }
}

mod generic_tests {
    use ark_bls12_381::Fr;
    use crate::data_structures::{GKRRound, Prover, Verifier};
    use crate::naive_sum_check::NaiveProver;
    use crate::standard_verifier::StandardVerifier;
    use crate::util::{random_gate};

    #[test]
    fn test_verifier_first_round() {
        // 7 variables 3 seconds.... 8 variables 23!!! seconds!??! :OOO whaaa

        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut prover = NaiveProver::new(gkr_round.clone(), &random_gate);
        // Now we test the g_func gives what we expect
        let verifier = StandardVerifier::new(3, prover.compute_sum(), gkr_round);
        let verifier_func = prover.get_verifier_function();

        assert!(verifier.check_claimed_value(&verifier_func));
    }

    #[test]
    fn simulate_two_rounds_naive() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let random_gate = random_gate(gkr_round.gate_labes());
        let mut prover = NaiveProver::new(gkr_round.clone(), &random_gate);
        // Now we test the g_func gives what we expect
        let mut verifier = StandardVerifier::new(3, prover.compute_sum(), gkr_round);
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