use std::time::{Duration, Instant};
use ark_bls12_381::Fr;
use ark_ff::Field;
use crate::data_structures::{GKRRound, SumCheckProver, SumCheckVerifier};
use crate::fast_prover::FastProver;
use crate::naive_sum_check::NaiveProver;
use crate::standard_verifier::StandardVerifier;
use crate::util::random_gate;

mod data_structures;
mod util;
pub mod naive_sum_check;
mod standard_verifier;
mod fast_prover;
pub mod circuit_structures;
pub mod gkr_protocol;

fn main() {
    simulate_two_rounds_fast();
}

fn simulate_two_rounds_fast() {
    let gkr_round: GKRRound<Fr> = GKRRound::new_rand_var_size(11);
    let random_gate = random_gate(gkr_round.gate_labes());
    let prover = FastProver::new(gkr_round.clone(), &random_gate);
    println!("The number of variables is: {:?}", gkr_round.gate_labes());
    println!("Running the test with fast prover");
    compare_verifier_sum(gkr_round.clone(), prover);
    println!("Running the test with naive prover");
    let prover = NaiveProver::new(gkr_round.clone(), &random_gate);
    compare_verifier_sum(gkr_round, prover);
}

fn compare_verifier_sum<F: Field, P: SumCheckProver<F>>(gkr_round: GKRRound<F>, mut prover: P) {
    let mut verifier = StandardVerifier::new(3, prover.compute_sum(), gkr_round.clone());
    let mut prover_time_spent = Duration::ZERO;
    let mut verifier_time_spent = Duration::ZERO;
    for _ in 0..gkr_round.gate_labes() {

        let time = Instant::now();
        let verifier_func = prover.get_verifier_function();
        let time_diff = time.elapsed();
        prover_time_spent += time_diff;

        let time = Instant::now();
        assert!(verifier.check_claimed_value(&verifier_func));
        let time_diff = time.elapsed();
        verifier_time_spent += time_diff;

        let time = Instant::now();
        let rand_var = verifier.get_random_field_element();
        let elapsed = time.elapsed();
        verifier_time_spent += elapsed;

        let time = Instant::now();
        prover.fix_variable(rand_var);
        let new_claim = prover.compute_sum();
        let time_diff = time.elapsed();
        prover_time_spent += time_diff;
        verifier.set_claim(new_claim);
    }
    println!("Prover time: {:?}", prover_time_spent);
    println!("Verifier time: {:?}", verifier_time_spent);
}

mod generic_tests {
    use std::time::{Duration, Instant};
    use ark_bls12_381::Fr;
    use ark_ff::Field;
    use ark_poly::MultilinearExtension;
    use crate::data_structures::{GKRRound, SumCheckProver, SumCheckVerifier};
    use crate::fast_prover::FastProver;
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
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand_var_size(10);
        let random_gate = random_gate(gkr_round.gate_labes());
        let prover = NaiveProver::new(gkr_round.clone(), &random_gate);
        compare_verifier_sum(gkr_round, prover);
    }

    #[test]
    fn simulate_two_rounds_fast() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand_var_size(8);
        let random_gate = random_gate(gkr_round.gate_labes());
        let prover = FastProver::new(gkr_round.clone(), &random_gate);
        compare_verifier_sum(gkr_round, prover);
        assert_eq!(1, 2);
    }

    fn compare_verifier_sum<F: Field, P: SumCheckProver<F>>(gkr_round: GKRRound<F>, mut prover: P) {
        let mut verifier = StandardVerifier::new(3, prover.compute_sum(), gkr_round.clone());
        let mut prover_time_spent = Duration::ZERO;
        let mut verifier_time_spent = Duration::ZERO;
        for _ in 0..gkr_round.gate_labes() {

            let verifier_func = prover.get_verifier_function();
            let time = Instant::now();
            assert!(verifier.check_claimed_value(&verifier_func));
            let time_diff = time.elapsed();
            verifier_time_spent += time_diff;

            let rand_var = verifier.get_random_field_element();
            prover.fix_variable(rand_var);

            let time = Instant::now();
            let new_claim = prover.compute_sum();
            let time_diff = time.elapsed();
            prover_time_spent += time_diff;
            verifier.set_claim(new_claim);
        }
        println!("Prover time: {:?}", prover_time_spent);
        println!("Verifier time: {:?}", verifier_time_spent);
    }

    #[test]
    fn test_generated_round() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand_var_size(8);
        assert_eq!(gkr_round.vi.num_vars, 8);
        assert_eq!(gkr_round.vj.num_vars, 8);
        assert_eq!(gkr_round.mult().num_vars(), 24);
    }
}