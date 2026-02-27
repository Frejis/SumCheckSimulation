use std::time::{Duration, Instant};
use ark_bls12_381::Fr;
use ark_ff::Field;
use structures::data_structures::{SumCheckProver, SumCheckVerifier};
use crate::gkr::gkr_round::GKRRound;
use crate::provers::{fast, naive};
use crate::verifiers::standard_verifier::StandardVerifier;
use crate::util::random_gate;

mod util;
pub mod provers;
pub mod structures;
pub mod verifiers;
pub mod gkr;

fn main() {
    simulate_and_print_results();
}

fn simulate_and_print_results() {
    for i in 0..10 {
        println!("##############################################################");
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand_var_size(i);
        let random_gate = random_gate(gkr_round.gate_labes());
        let prover = fast::FastProver::new(gkr_round.clone(), &random_gate);
        println!("The number of variables is: {:?}", gkr_round.gate_labes()*2);
        println!("Running the test with fast prover");
        compare_verifier_sum(gkr_round.clone(), prover);
        println!("Running the test with naive prover");
        let prover = naive::NaiveProver::new(gkr_round.clone(), &random_gate);
        compare_verifier_sum(gkr_round, prover);
    }
}


/// Prints the running time of the prover and verifier for a GKRRound.
///
/// # Arguments
///
/// * `gkr_round`:
/// * `prover`:
///
/// returns: ()
///
/// # Examples
///
/// ```
///     let gkr_round: GKRRound<Fr> = GKRRound::new_rand_var_size(11);
///     let random_gate = random_gate(gkr_round.gate_labes());
///     let prover = fast::FastProver::new(gkr_round.clone(), &random_gate);
///     println!("The number of variables is: {:?}", gkr_round.gate_labes());
///     println!("Running the test with fast prover");
///     compare_verifier_sum(gkr_round.clone(), prover);
///     println!("Running the test with naive prover");
///     let prover = naive::NaiveProver::new(gkr_round.clone(), &random_gate);
///     compare_verifier_sum(gkr_round, prover);
/// ```
fn compare_verifier_sum<F: Field, P: SumCheckProver<F>>(gkr_round: GKRRound<F>, mut prover: P) {
    let mut verifier = StandardVerifier::new(3, prover.compute_sum(), gkr_round.clone());
    let mut prover_time_spent = Duration::ZERO;
    let mut verifier_time_spent = Duration::ZERO;
    for _ in 0..gkr_round.gate_labes() {
        simulate_round_timed(&mut prover, &mut verifier, &mut prover_time_spent, &mut verifier_time_spent);
    }
    println!("Prover time: {:?}", prover_time_spent);
    println!("Verifier time: {:?}", verifier_time_spent);
}

fn simulate_round_timed<F: Field, P: SumCheckProver<F>>(prover: &mut P, verifier: &mut StandardVerifier<F>, prover_time_spent: &mut Duration, verifier_time_spent: &mut Duration) {
    let time = Instant::now();
    let verifier_func = prover.get_verifier_function();
    let time_diff = time.elapsed();
    *prover_time_spent += time_diff;

    let time = Instant::now();
    assert!(verifier.check_claimed_value(&verifier_func));
    let time_diff = time.elapsed();
    *verifier_time_spent += time_diff;

    let time = Instant::now();
    let rand_var = verifier.get_random_field_element();
    let elapsed = time.elapsed();
    *verifier_time_spent += elapsed;

    let time = Instant::now();
    prover.fix_variable(rand_var);
    let new_claim = prover.compute_sum();
    let time_diff = time.elapsed();
    *prover_time_spent += time_diff;
    verifier.set_claim(new_claim);
}

#[cfg(test)]
mod generic_tests {
    use std::time::{Duration, Instant};
    use ark_bls12_381::Fr;
    use ark_ff::Field;
    use ark_poly::MultilinearExtension;
    use crate::structures::data_structures::{SumCheckProver, SumCheckVerifier};
    use crate::fast::FastProver;
    use crate::gkr::gkr_round::GKRRound;
    use crate::naive::NaiveProver;
    use crate::util::random_gate;
    use crate::verifiers::standard_verifier::StandardVerifier;

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
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand_var_size(5);
        let random_gate = random_gate(gkr_round.gate_labes());
        let prover = NaiveProver::new(gkr_round.clone(), &random_gate);
        _compare_verifier_sum(gkr_round, prover);
    }

    #[test]
    fn simulate_two_rounds_fast() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand_var_size(5);
        let random_gate = random_gate(gkr_round.gate_labes());
        let prover = FastProver::new(gkr_round.clone(), &random_gate);
        _compare_verifier_sum(gkr_round, prover);
    }

    fn _compare_verifier_sum<F: Field, P: SumCheckProver<F>>(gkr_round: GKRRound<F>, mut prover: P) {
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