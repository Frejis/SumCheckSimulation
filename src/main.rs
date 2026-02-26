use std::time::{Duration, Instant};
use ark_bls12_381::Fr;
use ark_ff::Field;
use ark_std::test_rng;
use structures::data_structures::{GKRRound, SumCheckProver, SumCheckVerifier};
use crate::provers::{fast, naive};
use crate::provers::naive::NaiveProver;
use crate::structures::circuit_structures::{GKRCircuit, GateType};
use crate::structures::gkr_protocol::GKRLayerDriver;
use crate::verifiers::standard_verifier::StandardVerifier;
use crate::util::random_gate;

mod util;
pub mod provers;
pub mod structures;
pub mod verifiers;

fn main() {
    simulate_two_rounds_fast();
}

fn simulate_two_rounds_fast() {
    let circuit: GKRCircuit<Fr> = GKRCircuit::random(&[1, 2, 4, 8], &mut test_rng());

    let (naive_prove_time, verifier_naive_time) =
        GKRLayerDriver::simulate_gkr_circuit::<Fr, NaiveProver<Fr>, _>(
            circuit.clone(),
            |round_i, current_r| NaiveProver::new(round_i, &current_r),
        );

    let (fast_prove_time, verifier_fast_time) =
        GKRLayerDriver::simulate_gkr_circuit::<Fr, fast::FastProver<Fr>, _>(
            circuit,
            |round_i, current_r| fast::FastProver::new(round_i, &current_r),
        );

    println!("Naive Prover time: {:?}", naive_prove_time);
    println!("Fast Prover time: {:?}", fast_prove_time);
    println!("Naive Verifier time: {:?}", verifier_naive_time);
    println!("Fast Verifier time: {:?}", verifier_fast_time);
}

fn get_round_information(circuit: &GKRCircuit<Fr>, i: &usize) -> (usize, GKRRound<Fr>) {
    fn log2_pow2(n: usize) -> usize {
        assert!(n.is_power_of_two());
        n.trailing_zeros() as usize
    }

    let layer_i = &circuit.layers[*i];
    let layer_next = &circuit.layers[*i + 1];

    let k_i = log2_pow2(layer_i.values.len());
    let k_next = log2_pow2(layer_next.values.len());

    let (_add_i, mult_i) = layer_i.wiring_predicates(k_i, k_next);
    let w_next = layer_next.value_extension(k_next);

    let round_i = GKRRound::new(&mult_i, &w_next, &w_next, &GateType::Mul);
    (k_next, round_i)
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
    use crate::structures::data_structures::{GKRRound, SumCheckProver, SumCheckVerifier};
    use crate::fast::FastProver;
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
        compare_verifier_sum(gkr_round, prover);
    }

    #[test]
    fn simulate_two_rounds_fast() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand_var_size(5);
        let random_gate = random_gate(gkr_round.gate_labes());
        let prover = FastProver::new(gkr_round.clone(), &random_gate);
        compare_verifier_sum(gkr_round, prover);
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