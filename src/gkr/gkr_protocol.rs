use std::time::{Duration, Instant};
use ark_ff::Field;
use crate::gkr::gkr_round::GKRRound;
use crate::structures::circuit_structures::{GKRCircuit, GateType};
use crate::structures::data_structures::{SumCheckProver, SumCheckVerifier};
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

        for i in 0..rounds {
            let time_now = Instant::now();
            let g_j = prover.get_verifier_function();
            prover_time_spent += time_now.elapsed();

            let time_now = Instant::now();
            println!("Round {i}");
            let r_j = verifier.handle_round(&g_j);
            verifier_time_spent += time_now.elapsed();

            let time_now = Instant::now();
            prover.fix_variable(r_j);
            let new_claim = prover.compute_sum();
            prover_time_spent += time_now.elapsed();

            let time_now = Instant::now();
            verifier.set_claim(new_claim);
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

    // For each non-input layer i, prove claim about W_i(r_i) via layer i+1.
    for i in 0..(circuit.layers.len() - 1) {
        let layer_i = &circuit.layers[i];
        let layer_next = &circuit.layers[i + 1];

        let k_i = log2_pow2(layer_i.values.len());
        let k_next = log2_pow2(layer_next.values.len());

        let (add_i, mult_i) = layer_i.wiring_predicates(k_i, k_next);
        let w_next = layer_next.value_extension(k_next);

        let round_i = GKRRound::new(&mult_i, &add_i, &w_next, &w_next, &GateType::Mul);

        let time = Instant::now();
        let mut prover = prover_ctor(round_i.clone(), current_r.clone());
        prover_time_spent += time.elapsed();
        let time = Instant::now();
        let mut verifier = StandardVerifier::new(3, prover.compute_sum(), round_i);
        verifier_time_spent += time.elapsed();
        println!("Time taken to construct prover and verifier for layer {}: {:?}", i, time.elapsed());
        let ((next_r, _next_claim), p_time, v_time) =
            GKRLayerDriver::run_layer(&mut prover, &mut verifier, k_next);

        current_r = next_r;
        prover_time_spent += p_time;
        verifier_time_spent += v_time;
    }

    (prover_time_spent, verifier_time_spent)
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_poly::Polynomial;
    use ark_std::test_rng;
    use crate::gkr::gkr_protocol::GKRLayerDriver;
    use crate::gkr::gkr_round::GKRRound;
    use crate::provers::fast::FastProver;
    use crate::provers::naive::NaiveProver;
    use crate::structures::circuit_structures::GKRCircuit;
    use crate::structures::data_structures::SumCheckProver;
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

        let ((next_r, next_claim), _, _) = GKRLayerDriver::run_layer(&mut prover, &mut verifier, k);

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

            let (add_i, mult_i) = layer_i.wiring_predicates(k_i, k_next);
            let w_next = layer_next.value_extension(k_next);

            let round_i = GKRRound::new(&mult_i, &add_i, &w_next, &w_next, &GateType::Mul);

            let mut prover = NaiveProver::new(round_i.clone(), &current_r);
            let mut verifier = StandardVerifier::new(3, current_claim, round_i);
            println!("layer is: {i}");
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
    fn mult_circuit_simulation_two_naiveop() {
        use crate::structures::circuit_structures::GateType;

        fn log2_pow2(n: usize) -> usize {
            assert!(n.is_power_of_two());
            n.trailing_zeros() as usize
        }

        // mul-only for now (your circuit generator currently does this).
        let circuit: GKRCircuit<Fr> = GKRCircuit::random(&[1, 2, 4, 8, 16, 32, 64, 64*2], &mut test_rng());

        // FIRST: Verify the circuit is constructed correctly
        println!("\n=== Verifying Circuit Construction ===");
        for i in 0..(circuit.layers.len() - 1) {
            let layer_i = &circuit.layers[i];
            let layer_next = &circuit.layers[i + 1];

            println!("\nChecking layer {}: {} gates -> {} gates",
                     i, layer_i.values.len(), layer_next.values.len());

            // For each gate in layer i, verify its value matches the computation from layer i+1
            for (gate_idx, gate) in layer_i.gates.iter().enumerate() {
                let expected_val = match gate.typ {
                    GateType::Add => layer_next.values[gate.left] + layer_next.values[gate.right],
                    GateType::Mul => layer_next.values[gate.left] * layer_next.values[gate.right],
                };
                assert_eq!(layer_i.values[gate_idx], expected_val,
                           "Layer {} gate {} value mismatch", i, gate_idx);
            }
            println!("  ✓ All gates correctly computed");
        }

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

            println!("\n=== Processing Layer {} ===", i);
            println!("k_i={}, k_next={}, current_r.len()={}", k_i, k_next, current_r.len());
            println!("current_r = {:?}", current_r);
            println!("current_claim = {}", current_claim);

            let (add_i, mult_i) = layer_i.wiring_predicates(k_i, k_next);
            let w_next = layer_next.value_extension(k_next);

            // Debug: Check gate types
            let num_add = layer_i.gates.iter().filter(|g| matches!(g.typ, GateType::Add)).count();
            let num_mul = layer_i.gates.iter().filter(|g| matches!(g.typ, GateType::Mul)).count();
            println!("  Layer {} has {} Add gates, {} Mul gates", i, num_add, num_mul);

            let round_i = GKRRound::new(&mult_i, &add_i, &w_next, &w_next, &GateType::Mul);

            let mut prover = NaiveProver::new(round_i.clone(), &current_r);

            // The prover computes what sum_{b,c} mult(r,b,c)*W(b)*W(c) equals
            let sumcheck_computed = prover.compute_sum();
            println!("prover computed sum = {}", sumcheck_computed);

            // Also compute using the reference ark implementation
            let ark_sum = NaiveProver::ark_compute_sum_naive(&mult_i, &add_i, &w_next, &w_next, &current_r);
            println!("ark reference sum = {}", ark_sum);

            // Also compute directly what W_i(current_r) should be
            let w_i = layer_i.value_extension(k_i);
            let direct_eval = w_i.evaluate(&current_r);
            println!("direct W_i(r) eval = {}", direct_eval);

            // This MUST equal current_claim for the circuit to be sound
            if sumcheck_computed != current_claim {
                println!("ERROR: Mismatch at layer {}!", i);
                println!("  Expected (from claim): {}", current_claim);
                println!("  Got (from sumcheck):   {}", sumcheck_computed);
                println!("  Direct W_i(r):         {}", direct_eval);
                panic!("Circuit construction error: sum check doesn't match layer value");
            }

            // Create verifier with the CLAIM (what we're trying to prove)
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
            let layer_i = &circuit.layers[i];
            let layer_next = &circuit.layers[i + 1];

            let k_i = log2_pow2(layer_i.values.len());
            let k_next = log2_pow2(layer_next.values.len());

            let (add_i, mult_i) = layer_i.wiring_predicates(k_i, k_next);
            let w_next = layer_next.value_extension(k_next);

            let round_i = GKRRound::new(&mult_i, &add_i, &w_next, &w_next, &GateType::Mul);

            let mut prover = FastProver::new(round_i.clone(), &current_r);
            let mut verifier = StandardVerifier::new(3, current_claim, round_i);
            println!("layer is: {i}");
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
    fn mult_circuit_simulation_two_fastop() {
        use crate::structures::circuit_structures::GateType;

        fn log2_pow2(n: usize) -> usize {
            assert!(n.is_power_of_two());
            n.trailing_zeros() as usize
        }

        // mul-only for now (your circuit generator currently does this).
        let circuit: GKRCircuit<Fr> = GKRCircuit::random(&[1, 2, 4, 8, 16, 32, 64, 64*2], &mut test_rng());

        // FIRST: Verify the circuit is constructed correctly
        println!("\n=== Verifying Circuit Construction ===");
        for i in 0..(circuit.layers.len() - 1) {
            let layer_i = &circuit.layers[i];
            let layer_next = &circuit.layers[i + 1];

            println!("\nChecking layer {}: {} gates -> {} gates",
                     i, layer_i.values.len(), layer_next.values.len());

            // For each gate in layer i, verify its value matches the computation from layer i+1
            for (gate_idx, gate) in layer_i.gates.iter().enumerate() {
                let expected_val = match gate.typ {
                    GateType::Add => layer_next.values[gate.left] + layer_next.values[gate.right],
                    GateType::Mul => layer_next.values[gate.left] * layer_next.values[gate.right],
                };
                assert_eq!(layer_i.values[gate_idx], expected_val,
                           "Layer {} gate {} value mismatch", i, gate_idx);
            }
            println!("  ✓ All gates correctly computed");
        }

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

            println!("\n=== Processing Layer {} ===", i);
            println!("k_i={}, k_next={}, current_r.len()={}", k_i, k_next, current_r.len());
            println!("current_r = {:?}", current_r);
            println!("current_claim = {}", current_claim);

            let (add_i, mult_i) = layer_i.wiring_predicates(k_i, k_next);
            let w_next = layer_next.value_extension(k_next);

            // Debug: Check gate types
            let num_add = layer_i.gates.iter().filter(|g| matches!(g.typ, GateType::Add)).count();
            let num_mul = layer_i.gates.iter().filter(|g| matches!(g.typ, GateType::Mul)).count();
            println!("  Layer {} has {} Add gates, {} Mul gates", i, num_add, num_mul);

            let round_i = GKRRound::new(&mult_i, &add_i, &w_next, &w_next, &GateType::Mul);

            let mut prover = FastProver::new(round_i.clone(), &current_r);

            // The prover computes what sum_{b,c} mult(r,b,c)*W(b)*W(c) equals
            let sumcheck_computed = prover.compute_sum();
            println!("prover computed sum = {}", sumcheck_computed);

            // Also compute using the reference ark implementation
            let ark_sum = NaiveProver::ark_compute_sum_naive(&mult_i, &add_i, &w_next, &w_next, &current_r);
            println!("ark reference sum = {}", ark_sum);

            // Also compute directly what W_i(current_r) should be
            let w_i = layer_i.value_extension(k_i);
            let direct_eval = w_i.evaluate(&current_r);
            println!("direct W_i(r) eval = {}", direct_eval);

            // This MUST equal current_claim for the circuit to be sound
            if sumcheck_computed != current_claim {
                println!("ERROR: Mismatch at layer {}!", i);
                println!("  Expected (from claim): {}", current_claim);
                println!("  Got (from sumcheck):   {}", sumcheck_computed);
                println!("  Direct W_i(r):         {}", direct_eval);
                panic!("Circuit construction error: sum check doesn't match layer value");
            }

            // Create verifier with the CLAIM (what we're trying to prove)
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
}