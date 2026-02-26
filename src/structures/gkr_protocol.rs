use std::time::{Duration, Instant};
use ark_bls12_381::Fr;
use ark_ff::Field;
use ark_poly::Polynomial;
use crate::provers::naive::NaiveProver;
use crate::structures::circuit_structures::{GateType, GKRCircuit};
use crate::structures::data_structures::{GKRRound, SumCheckProver, SumCheckVerifier};
use crate::util::random_gate;
use crate::verifiers::standard_verifier::StandardVerifier;

pub struct GKRLayerDriver;

impl GKRLayerDriver {
    pub fn run_layer<F, P> (
        prover: &mut P,
        verifier: &mut StandardVerifier<F>,
        k_next: usize,
    ) -> ((Vec<F>, F), Duration, Duration)
    where
        F: Field,
        P: SumCheckProver<F>,
    {
        let rounds = 2 * k_next;

        let mut prover_time_spent = Duration::ZERO;
        let mut verifier_time_spent = Duration::ZERO;

        for _ in 0..rounds {
            let time_now = Instant::now();
            let g_j = prover.get_verifier_function();
            prover_time_spent += time_now.elapsed();
            let time_now = Instant::now();
            let r_j = verifier.handle_round(&g_j);
            verifier_time_spent += time_now.elapsed();
            let time_now = Instant::now();
            prover.fix_variable(r_j);
            verifier_time_spent += time_now.elapsed();

            let time_now = Instant::now();
            let next_claim_for_round = prover.compute_sum();
            verifier_time_spent += time_now.elapsed();

            let time_now = Instant::now();
            verifier.set_claim(next_claim_for_round);
            verifier_time_spent += time_now.elapsed();
        }

        let points = verifier.random_points_chosen();
        assert!(points.len() >= rounds, "not enough sampled points in verifier transcript");

        let tail = &points[points.len() - rounds..]; // Note to self this is because the verifier will always push
        // So points will at some point include ALL variables chosen by the verifier.

        let time_now = Instant::now();
        let b_star = tail[0..k_next].to_vec();
        let c_star = tail[k_next..2 * k_next].to_vec();

        assert_eq!(b_star.len(), k_next, "b_star should have length k_next");
        assert_eq!(c_star.len(), k_next, "c_star should have length k_next");

        let msg = prover.layer_reduction_message(&b_star, &c_star);
        prover_time_spent += time_now.elapsed();

        let mut rng = ark_std::test_rng();
        let time_now = Instant::now();
        let left_res = verifier.reduce_two_claims_to_one(&b_star, &c_star, &msg, &mut rng);
        verifier_time_spent += time_now.elapsed();
        (left_res, prover_time_spent, verifier_time_spent)
    }

    pub fn simulate_gkr_circuit<F: Field, P, C>(
        circuit: GKRCircuit<F>,
        mut prover_ctor: C,
    ) -> (Duration, Duration)
    where
        P: SumCheckProver<F>,
        C: FnMut(GKRRound<F>, Vec<F>) -> P,
    {
        fn log2_pow2(n: usize) -> usize {
            assert!(n.is_power_of_two());
            n.trailing_zeros() as usize
        }

        // Time prover is the first in the result and the other is verifier.
        let mut prover_time_spent = Duration::ZERO;
        let mut verifier_time_spent = Duration::ZERO;

        let k0 = log2_pow2(circuit.layers[0].values.len());

        let mut current_r = random_gate::<F>(k0);
        let mut current_claim = circuit.layers[0].value_extension(k0).evaluate(&current_r);

        // For each non-input layer i, prove claim about W_i(r_i) via layer i+1.
        for i in 0..(circuit.layers.len() - 1) {
            let layer_i = &circuit.layers[i];
            let layer_next = &circuit.layers[i + 1];

            let k_i = log2_pow2(layer_i.values.len());
            let k_next = log2_pow2(layer_next.values.len());

            let (_add_i, mult_i) = layer_i.wiring_predicates(k_i, k_next);
            let w_next = layer_next.value_extension(k_next);

            let round_i = GKRRound::new(&mult_i, &w_next, &w_next, &GateType::Mul);

            let mut prover = prover_ctor(round_i.clone(), current_r.clone());
            let mut verifier = StandardVerifier::new(3, current_claim, round_i);

            let ((next_r, next_claim), p_time, v_time) =
                GKRLayerDriver::run_layer(&mut prover, &mut verifier, k_next);

            current_r = next_r;
            current_claim = next_claim;
            prover_time_spent += p_time;
            verifier_time_spent += v_time;
        }

        (prover_time_spent, verifier_time_spent)
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_poly::{MultilinearExtension, Polynomial};
    use ark_std::test_rng;
    use crate::provers::fast::FastProver;
    use crate::structures::circuit_structures::{GateType, GKRCircuit};
    use crate::structures::data_structures::{GKRRound, SumCheckProver};
    use crate::provers::naive::NaiveProver;
    use crate::structures::gkr_protocol::GKRLayerDriver;
    use crate::util::random_gate;
    use crate::verifiers::standard_verifier::StandardVerifier;

    #[test]
    fn simulate_full_gkr_round_with_naive_prover() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let k = gkr_round.vi().num_vars;
        let fixed_gate = random_gate::<Fr>(gkr_round.gate_labes());

        let mut prover = NaiveProver::new(gkr_round.clone(), &fixed_gate);

        // IMPORTANT: verifier initial claim must match prover's current sum.
        let initial_claim = prover.compute_sum();
        let mut verifier = StandardVerifier::new(3, initial_claim, gkr_round.clone());

        let ((next_r, next_claim), _, _)  = GKRLayerDriver::run_layer(&mut prover, &mut verifier, k);

        // Folded claim must match direct W_{i+1}(r_next) evaluation.
        let expected = gkr_round.vi().evaluate(&next_r);
        assert_eq!(next_claim, expected);
    }

    #[test]
    fn mult_circuit_simulation_two_naive() {
        use crate::structures::circuit_structures::GateType;

        fn log2_pow2(n: usize) -> usize {
            assert!(n.is_power_of_two());
            n.trailing_zeros() as usize
        }

        // mul-only for now (your circuit generator currently does this).
        let circuit: GKRCircuit<Fr> = GKRCircuit::random(&[1, 2, 4, 8, 16, 32, 64, 64*2], &mut test_rng());

        // Start with random output-layer point r_0.
        let k0 = log2_pow2(circuit.layers[0].values.len());
        let mut current_r = random_gate::<Fr>(k0);
        let mut current_claim = circuit.layers[0].value_extension(k0).evaluate(&current_r);

        // For each non-input layer i, prove claim about W_i(r_i) via layer i+1.
        for i in 0..(circuit.layers.len() - 1) {
            let layer_i = &circuit.layers[i];
            let layer_next = &circuit.layers[i + 1];

            let k_i = log2_pow2(layer_i.values.len());
            let k_next = log2_pow2(layer_next.values.len());

            let (_add_i, mult_i) = layer_i.wiring_predicates(k_i, k_next);
            let w_next = layer_next.value_extension(k_next);

            let round_i = GKRRound::new(&mult_i, &w_next, &w_next, &GateType::Mul);

            let mut prover = NaiveProver::new(round_i.clone(), &current_r);
            let mut verifier = StandardVerifier::new(3, current_claim, round_i);

            let ((next_r, next_claim), _, _)  = GKRLayerDriver::run_layer(&mut prover, &mut verifier, k_next);

            current_r = next_r;
            current_claim = next_claim;
        }

        // Final check at input layer.
        let last_idx = circuit.layers.len() - 1;
        let k_last = log2_pow2(circuit.layers[last_idx].values.len());
        let w_last = circuit.layers[last_idx].value_extension(k_last);
        let expected = w_last.evaluate(&current_r);

        assert_eq!(current_claim, expected, "Final input-layer claim mismatch");
    }

    #[test]
    fn mult_circuit_simulation_two_fast() {
        use crate::structures::circuit_structures::GateType;

        fn log2_pow2(n: usize) -> usize {
            assert!(n.is_power_of_two());
            n.trailing_zeros() as usize
        }

        // mul-only for now (your circuit generator currently does this).
        let circuit: GKRCircuit<Fr> = GKRCircuit::random(&[1, 2, 4, 8, 16, 32, 64, 64*2], &mut test_rng());

        // Start with random output-layer point r_0.
        let k0 = log2_pow2(circuit.layers[0].values.len());
        let mut current_r = random_gate::<Fr>(k0);
        let mut current_claim = circuit.layers[0].value_extension(k0).evaluate(&current_r);

        // For each non-input layer i, prove claim about W_i(r_i) via layer i+1.
        for i in 0..(circuit.layers.len() - 1) {
            let (k_next, round_i) = get_round_information(&circuit, &i);

            let mut prover = FastProver::new(round_i.clone(), &current_r);
            let mut verifier = StandardVerifier::new(3, current_claim, round_i);

            let ((next_r, next_claim), _, _)  = GKRLayerDriver::run_layer(&mut prover, &mut verifier, k_next);

            current_r = next_r;
            current_claim = next_claim;
        }

        // Final check at input layer.
        let last_idx = circuit.layers.len() - 1;
        let k_last = log2_pow2(circuit.layers[last_idx].values.len());
        let w_last = circuit.layers[last_idx].value_extension(k_last);
        let expected = w_last.evaluate(&current_r);

        assert_eq!(current_claim, expected, "Final input-layer claim mismatch");
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
}