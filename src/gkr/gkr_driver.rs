use std::time::{Duration};
use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, Polynomial};
use crate::gkr::gkr_prover::GKRProver;
use crate::gkr::gkr_round::GKRRound;
use crate::gkr::gkr_verifier::GKRVerifier;
use crate::gkr::layer::InputLayer;
use crate::provers::fast::FastProver;
use crate::structures::circuit_structures::{GKRCircuit};
use crate::structures::data_structures::{SumCheckProver, SumCheckVerifier};
use crate::verifiers::standard_verifier::StandardVerifier;

/// This file is responsible for simulating a GKR Proof.
/// It implements a driver that takes a Circuit, Prover, Verifier and can simulate the entire 
/// GKR protocol. There are two versions of this. 
/// Currently, it does not benchmark the time prover/verfier spent.

pub struct GKRDriver<F: Field> {
    gkrprover: GKRProver<F>,
    verifier: GKRVerifier<F>,
    circuit: GKRCircuit<F>,
    input_layer: InputLayer<F>,
}

impl<F: Field> GKRDriver<F> {
    pub fn new(gkrprover: GKRProver<F>, verifier: GKRVerifier<F>, circuit: GKRCircuit<F>, input_layer: InputLayer<F>) -> Self {
        Self { gkrprover, verifier, circuit, input_layer }
    }
}

pub fn log2_pow2(n: usize) -> usize {
    assert!(n.is_power_of_two());
    n.trailing_zeros() as usize
}



impl<F: Field> GKRDriver<F> {

    /// Runs sum check for s_i_plus_1 rounds.
    /// Returns a random gate chosen via reducing to claims to a claim about one
    /// Returns the alleged claim for that position.
    pub fn run_layer (
        mut prover: FastProver<F>,
        mut verifier: StandardVerifier<F>,
        s_i_plus_1: usize,
    ) -> (Vec<F>, F)
    {
        let total_rounds = 2 * s_i_plus_1;
        for _ in 0..total_rounds {
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
    ) -> (Duration, Duration)
    {
        let mut mi = F::zero();
        let mut next_gate = vec![F::zero()];
        // This is for all rounds except the last round.
        for i in 0..self.circuit.layers.len() {
            let s_i_plus_1 = self.get_correct_next_layer_size(i);
            let value_extension = self.get_correct_value_extension(i);
            if i == 0 {
                println!("Running initial round");
                (next_gate, mi) = self.handle_first_round(&mut next_gate, s_i_plus_1, &value_extension);
            } else {
                println!("Handling other rounds");
                (next_gate, mi) = self.handle_intermediate_rounds(mi, &mut next_gate, s_i_plus_1, &value_extension, i);
            }
        }

        // This function panics == Verifier rejects.
        println!("Checking Final claim");
        self.verifier.verify_final_claimed_value_point(next_gate, mi);

        (Duration::from_hours(2), Duration::from_hours(2))
    }

    fn get_correct_next_layer_size(&mut self, i: usize) -> usize {
        if i < self.circuit.layers.len() - 1 {
            let layer = &self.circuit.layers[i + 1];
            let gates_len = layer.gates.len();
            log2_pow2(gates_len)
        } else {
            log2_pow2(self.input_layer.values.len())
        }
    }

    fn get_correct_value_extension(&mut self, i: usize) -> DenseMultilinearExtension<F> {
        if i < self.circuit.layers.len() - 1 {
            self.gkrprover
                .evaluated_circuit()
                .layers[i + 1]
                .value_extension()
        } else {
            self.input_layer
                .value_extension()
        }
    }

    fn handle_intermediate_rounds(&mut self,
                                  mi: F,
                                  next_gate: &mut Vec<F>,
                                  s_i_plus_1: usize,
                                  value_extension: &DenseMultilinearExtension<F>,
                                  layer: usize,
    ) -> (Vec<F>, F) {
        let (add_pred, mult_pred) = &self.gkrprover.predicates()[layer];
        let gkr_round: GKRRound<F> = GKRRound::new(&mult_pred.pred, &add_pred.pred, &value_extension, &value_extension);
        let mut fast_prover = FastProver::new(gkr_round.clone(), &*next_gate);
        assert_eq!(mi, fast_prover.compute_sum());
        let verifier = StandardVerifier::new(100, mi, gkr_round);
        Self::run_layer(fast_prover, verifier, s_i_plus_1)
    }

    fn handle_first_round(&mut self,
                          next_gate: &mut Vec<F>,
                          s_i_plus_1: usize,
                          value_extension: &DenseMultilinearExtension<F>
    ) -> (Vec<F>, F) {
        let output_claim = self.gkrprover.get_output_claim();
        let (add_pred, mult_pred) = &self.gkrprover.predicates()[0];
        let gkr_round: GKRRound<F> = GKRRound::new(&mult_pred.pred, &add_pred.pred, &value_extension.clone(), &value_extension.clone());
        let mut fast_prover = FastProver::new(gkr_round.clone(), &*next_gate);

        assert_eq!(output_claim.evaluate(&next_gate), fast_prover.compute_sum()); // This will be kept here as a sanity check.


        let mut fast_prover = FastProver::new(gkr_round.clone(), &*next_gate);
        assert_eq!(output_claim.evaluate(&next_gate), fast_prover.compute_sum());
        let m0 = fast_prover.compute_sum();
        let verifier = StandardVerifier::new(100, m0, gkr_round);
        Self::run_layer(fast_prover, verifier, s_i_plus_1)
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_poly::Polynomial;
    use crate::gkr::gkr_driver::GKRDriver;
    use crate::gkr::gkr_round::GKRRound;
    use crate::provers::fast::FastProver;
    use crate::structures::data_structures::SumCheckProver;
    use crate::util::random_gate;
    use crate::verifiers::standard_verifier::StandardVerifier;

    #[test]
    fn simulate_full_gkr_round_with_naive_prover() {
        let gkr_round: GKRRound<Fr> = GKRRound::new_rand(7);
        let k = gkr_round.vi().num_vars;
        let fixed_gate = random_gate::<Fr>(gkr_round.gate_labes());

        let mut prover = FastProver::new(gkr_round.clone(), &fixed_gate);

        // IMPORTANT: verifier initial claim must match prover's current sum.
        let initial_claim = prover.compute_sum();
        let verifier = StandardVerifier::new(3, initial_claim, gkr_round.clone());

        let (next_r, next_claim) = GKRDriver::run_layer(prover, verifier, k);

        // Folded claim must match direct W_{i+1}(r_next) evaluation.
        let expected = gkr_round.vi().evaluate(&next_r);
        assert_eq!(next_claim, expected);
    }


}