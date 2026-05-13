use std::time::{Duration, Instant};
use ark_ff::Field;
use crate::gkr::gkr_prover::GKRProver;
use crate::gkr::gkr_round::GKRRound;
use crate::gkr::gkr_verifier::GKRVerifier;
use crate::provers::fast::FastProver;
use crate::structures::circuit_structures::GKRCircuit;
use crate::structures::data_structures::{SumCheckProver, SumCheckVerifier};
use crate::verifiers::standard_verifier::StandardVerifier;

/// This file is responsible for simulating a GKRProof
/// It implements a driver that takes a Circuit, Prover, Verifier and can simulate the entire 
/// GKR protocol. There are two versions of this. 
/// Currently it does not benchmark the time prover/verfier spent.

struct GKRDriver<F: Field> {
    GKRProver: GKRProver<F>,
    verifier: GKRVerifier<F>,
    circuit: GKRCircuit<F>,
}

pub fn log2_pow2(n: usize) -> usize {
    assert!(n.is_power_of_two());
    n.trailing_zeros() as usize
}



impl<F: Field> GKRDriver<F> {
    pub fn run_layer (
        mut prover: FastProver<F>,
        mut verifier: StandardVerifier<F>,
        s_i_plus_1: usize,
    ) -> (Vec<F>, F)
    {
        let total_rounds = 2 * s_i_plus_1;
        for i in 0..total_rounds {
            let g_j = prover.get_verifier_function();

            let r_j = verifier.handle_round(&g_j);

            prover.fix_variable(r_j);

            let new_claim = prover.compute_sum();

            verifier.set_claim(new_claim);
        }

        let msg = prover.layer_reduction_message(s_i_plus_1);
        let (next_layer_gate, claim) = verifier.handle_layer_reduction_message(msg, s_i_plus_1);
        (next_layer_gate, claim)
    }

    pub fn run_circuit (
        &mut self,
        k_next: usize,
    ) -> (Duration, Duration)
    {
        let mut fixed_gate: &[F] = &vec![];
        let mut mi = F::zero();
        let mut next_gate = vec![F::zero()];
        // This is for all rounds except the last round.
        for i in 0..self.circuit.layers.len() - 1 {
            let s_i_plus_1 = log2_pow2(self.circuit.layers[i+1].gates.len());
            let value_extension = self.GKRProver
                .evaluated_circuit()
                .layers[i]
                .value_extension(s_i_plus_1);

            if i == 0 {
                // TODO: Handle round 1 where we get the initial claimed output and fix the first gate.
                let output_claim = self.GKRProver.get_output_claim();
                let gate_length = log2_pow2(self.circuit.layers[0].gates.len());
                let random_gate = self.verifier.random_gate(&output_claim, gate_length);
                let (add_pred, mult_pred) = &self.GKRProver.predicates()[0];
                let gkr_round: GKRRound<F> = GKRRound::new(&add_pred.pred, &mult_pred.pred, &value_extension, &value_extension);
                let mut fast_prover = FastProver::new(gkr_round.clone(), &*random_gate);
                let m0 = fast_prover.compute_sum();
                let verifier = StandardVerifier::new(100, m0, gkr_round);
                (next_gate, mi) = Self::run_layer(fast_prover, verifier, s_i_plus_1);
            } else {
                let (add_pred, mult_pred) = &self.GKRProver.predicates()[0];
                let gkr_round: GKRRound<F> = GKRRound::new(&add_pred.pred, &mult_pred.pred, &value_extension, &value_extension);
                let mut fast_prover = FastProver::new(gkr_round.clone(), &*next_gate);
                let verifier = StandardVerifier::new(100, mi, gkr_round);
                (next_gate, mi) = Self::run_layer(fast_prover, verifier, s_i_plus_1);
            }
        }
        (Duration::from_hours(2), Duration::from_hours(2))
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_poly::Polynomial;
    use crate::gkr::gkr_driver::GKRDriver;
    use crate::gkr::gkr_round::GKRRound;
    use crate::provers::fast::FastProver;
    use crate::provers::naive::NaiveProver;
    use crate::structures::data_structures::SumCheckProver;
    use crate::util::random_gate;
    use crate::verifiers::standard_verifier::StandardVerifier;

    #[test]
    fn simulate_full_gkr_round_with_naive_prover() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand();
        let k = gkr_round.vi().num_vars;
        let fixed_gate = random_gate::<Fr>(gkr_round.gate_labes());

        let mut prover = FastProver::new(gkr_round.clone(), &fixed_gate);

        // IMPORTANT: verifier initial claim must match prover's current sum.
        let initial_claim = prover.compute_sum();
        let mut verifier = StandardVerifier::new(3, initial_claim, gkr_round.clone());

        let ((next_r, next_claim)) = GKRDriver::run_layer(prover, verifier, k);

        // Folded claim must match direct W_{i+1}(r_next) evaluation.
        let expected = gkr_round.vi().evaluate(&next_r);
        assert_eq!(next_claim, expected);
    }


}