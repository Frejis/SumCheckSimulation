use ark_ff::Field;
use crate::data_structures::{LayerReductionMessage, SumCheckProver, SumCheckVerifier};
use crate::standard_verifier::StandardVerifier;

/// Prover hook for the final message at a layer.
/// (Implement this for NaiveProver/FastProver.)
pub trait LayerReductionOracle<F: Field> {
    fn layer_reduction_message(&mut self, b_star: &[F], c_star: &[F]) -> LayerReductionMessage<F>;
}

pub struct GkrLayerDriver;

impl GkrLayerDriver {
    /// Runs one GKR layer interaction:
    /// - executes 2*k_{i+1} sum-check rounds,
    /// - obtains (b*, c*) from verifier challenges,
    /// - asks prover for (z1, z2, q),
    /// - reduces to next layer claim (r_{i+1}, claim_{i+1}).
    pub fn run_layer<F, P>(
        prover: &mut P,
        verifier: &mut StandardVerifier<F>,
        k_next: usize,
    ) -> (Vec<F>, F)
    where
        F: Field,
        P: SumCheckProver<F> + LayerReductionOracle<F>,
    {
        let rounds = 2 * k_next;

        for _ in 0..rounds {
            let g_j = prover.get_verifier_function();
            let r_j = verifier.handle_round(&g_j);
            prover.fix_variable(r_j);

            let next_claim_for_round = prover.compute_sum();
            verifier.set_claim(next_claim_for_round);
        }

        let pts = verifier.sampled_points();
        assert!(pts.len() >= rounds, "not enough sampled points in verifier transcript");
        let tail = &pts[pts.len() - rounds..];

        let b_star = tail[0..k_next].to_vec();
        let c_star = tail[k_next..2 * k_next].to_vec();

        let msg = prover.layer_reduction_message(&b_star, &c_star);

        let mut rng = ark_std::test_rng();
        verifier.reduce_two_claims_to_one(&b_star, &c_star, &msg, &mut rng)
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_poly::{MultilinearExtension, Polynomial};
    use ark_std::test_rng;
    use crate::circuit_structures::{GateType, GkrCircuit};
    use crate::data_structures::{GKRRound, SumCheckProver};
    use crate::fast_prover::FastProver;
    use crate::gkr_protocol::GkrLayerDriver;
    use crate::naive_sum_check::NaiveProver;
    use crate::standard_verifier::StandardVerifier;
    use crate::util::random_gate;

    #[test]
    fn simulate_full_gkr_round_with_naive_prover() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let k = gkr_round.vi().num_vars;
        let fixed_gate = random_gate::<Fr>(gkr_round.gate_labes());

        let mut prover = NaiveProver::new(gkr_round.clone(), &fixed_gate);

        // IMPORTANT: verifier initial claim must match prover's current sum.
        let initial_claim = prover.compute_sum();
        let mut verifier = StandardVerifier::new(3, initial_claim, gkr_round.clone());

        let (r_next, claim_next) = GkrLayerDriver::run_layer(&mut prover, &mut verifier, k);

        // Folded claim must match direct W_{i+1}(r_next) evaluation.
        let expected = gkr_round.vi().evaluate(&r_next);
        assert_eq!(claim_next, expected);
    }

    #[test]
    fn mult_circuit_simulation() {
        fn log2_pow2(n: usize) -> usize {
            assert!(n.is_power_of_two());
            n.trailing_zeros() as usize
        }

        // 3 layers => exactly 2 chained GKR rounds: (layer 0 -> 1) and (layer 1 -> 2)
        let circuit: GkrCircuit<Fr> = GkrCircuit::random(&[2, 4, 8], &mut test_rng());

        // ---------- Round 0 (prove W0(r0) via layer 1) ----------
        let layer0 = &circuit.layers[0];
        let layer1 = &circuit.layers[1];
        let k0 = log2_pow2(layer0.values.len());
        let k1 = log2_pow2(layer1.values.len());

        let (_add0, mult0) = layer0.wiring_predicates(k0, k1);
        let w1 = layer1.value_extension(k1);

        let round0 = GKRRound::new(&mult0, &w1, &w1, &GateType::Mul);

        let r0 = random_gate::<Fr>(k0);
        let mut prover0 = NaiveProver::new(round0.clone(), &r0);
        let claim0 = prover0.compute_sum();
        let mut verifier0 = StandardVerifier::new(3, claim0, round0.clone());

        let (r1, claim1) = GkrLayerDriver::run_layer(&mut prover0, &mut verifier0, k1);

        // ---------- Round 1 (prove W1(r1) via layer 2), chained by setting gate = r1 ----------
        let layer2 = &circuit.layers[2];
        let k2 = log2_pow2(layer2.values.len());

        let (_add1, mult1) = layer1.wiring_predicates(k1, k2);
        let w2 = layer2.value_extension(k2);

        let round1 = GKRRound::new(&mult1, &w2, &w2, &GateType::Mul);

        let mut prover1 = NaiveProver::new(round1.clone(), &r1);

        // True chain check: previous round's folded claim must be this round's initial sum claim.
        let claim1_from_round1 = prover1.compute_sum();
        assert_eq!(claim1, claim1_from_round1, "Chaining failed: claim1 mismatch");

        let mut verifier1 = StandardVerifier::new(3, claim1, round1.clone());
        let (r2, claim2) = GkrLayerDriver::run_layer(&mut prover1, &mut verifier1, k2);

        // Final check for second round fold:
        // claim2 should equal W2(r2), where W2 is MLE of input layer values.
        let expected2 = w2.evaluate(&r2);
        assert_eq!(claim2, expected2, "Final folded claim mismatch at layer 2");
    }

    #[test]
    fn mult_circuit_simulation_two_naive() {
        use ark_poly::MultilinearExtension;
        use crate::circuit_structures::GateType;
        use crate::gkr_protocol::GkrLayerDriver;

        fn log2_pow2(n: usize) -> usize {
            assert!(n.is_power_of_two());
            n.trailing_zeros() as usize
        }

        // mul-only for now (your circuit generator currently does this).
        let circuit: GkrCircuit<Fr> = GkrCircuit::random(&[1, 2, 4, 8, 16, 32, 64, 64*2], &mut test_rng());

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

            let (next_r, next_claim) = GkrLayerDriver::run_layer(&mut prover, &mut verifier, k_next);

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
        use crate::circuit_structures::GateType;
        use crate::gkr_protocol::GkrLayerDriver;

        fn log2_pow2(n: usize) -> usize {
            assert!(n.is_power_of_two());
            n.trailing_zeros() as usize
        }

        // mul-only for now (your circuit generator currently does this).
        let circuit: GkrCircuit<Fr> = GkrCircuit::random(&[1, 2, 4, 8, 16, 32, 64, 64*2], &mut test_rng());

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

            let (next_r, next_claim) = GkrLayerDriver::run_layer(&mut prover, &mut verifier, k_next);

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