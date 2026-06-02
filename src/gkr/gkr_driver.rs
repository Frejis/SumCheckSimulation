use std::time::{Duration, Instant};
use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, Polynomial, SparseMultilinearExtension};
use ark_poly::univariate::SparsePolynomial;
use crate::gkr::gkr_prover::GKRProver;
use crate::gkr::gkr_round::GKRRound;
use crate::gkr::gkr_verifier::GKRVerifier;
use crate::gkr::layer::{InputLayer, LayerConnection, LayerReductionMessage};
use crate::structures::circuit_structures::{GKRCircuit};
use crate::structures::data_structures::{AnalysisResult, SumCheckProver, SumCheckVerifier, Track};
use crate::util::create_prover;
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
    pub fn new(gkrprover: &GKRProver<F>, verifier: &GKRVerifier<F>, circuit: GKRCircuit<F>, input_layer: InputLayer<F>) -> Self {
        Self { gkrprover: gkrprover.clone(), verifier: verifier.clone(), circuit, input_layer }
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
    pub fn run_layer<T: SumCheckProver<F>> (
        &mut self,
        mut prover: T,
        mut verifier: StandardVerifier<F>,
        s_i_plus_1: usize,
        layer: usize,
    ) -> (LayerConnection<F>, Track)
    {
        let mut prover_time = Duration::new(0, 0);
        let mut verifier_time = Duration::new(0,0);
        let mut points;
        let total_rounds = 2 * s_i_plus_1;
        println!("Total rounds {total_rounds}");
        for _ in 0..total_rounds {
            points = Self::track_get_verif_func(&mut prover, &mut prover_time);

            // r_j is the random element returned by verf. in round j of sum-check.
            let r_j = Self::handle_verifier_function_from_prover(&mut verifier, &mut verifier_time, points);

            // Prover
            let new_claim = Self::track_fix_variable_and_new_claim(&mut prover, &mut prover_time, r_j);

            // Verifier
            Self::track_setting_new_claim(&mut verifier, &mut verifier_time, new_claim);
        }
        let msg = Self::track_creating_layer_reduc_msg(&mut prover, s_i_plus_1, &mut prover_time);

        let (next_layer_gate, claim) = self.track_handle_layer_reduc_msg(&mut verifier, s_i_plus_1, &mut verifier_time, msg, layer);

        self.verifier.set_next_layer_claim(claim);

        (LayerConnection::new(next_layer_gate, claim), Track::new_times(prover_time, verifier_time))
    }

    fn track_get_verif_func<T: SumCheckProver<F>>(prover: &mut T, prover_time: &mut Duration) -> Vec<(usize, F)> {
        let inst = Instant::now();
        let points = prover.get_verifier_function();
        let elapsed = inst.elapsed();
        *prover_time += elapsed;
        points
    }

    fn handle_verifier_function_from_prover(verifier: &mut StandardVerifier<F>, verifier_time: &mut Duration, points: Vec<(usize, F)>) -> F {
        let inst = Instant::now();
        let g_j = SparsePolynomial::from_coefficients_vec(points);
        let r_j = verifier.handle_round(&g_j);
        *verifier_time += inst.elapsed();
        r_j
    }

    fn track_setting_new_claim(verifier: &mut StandardVerifier<F>, verifier_time: &mut Duration, new_claim: F) {
        let inst = Instant::now();
        verifier.set_claim(new_claim);
        *verifier_time += inst.elapsed();
    }

    fn track_fix_variable_and_new_claim<T: SumCheckProver<F>>(prover: &mut T, prover_time: &mut Duration, r_j: F) -> F {
        let inst = Instant::now();
        prover.fix_variable(r_j);
        let new_claim = prover.compute_sum();
        let elapsed = inst.elapsed();
        *prover_time += elapsed;
        println!("Time taken to compute sum: {:?}", elapsed);
        new_claim
    }

    fn track_handle_layer_reduc_msg(&mut self, verifier: &mut StandardVerifier<F>, s_i_plus_1: usize, verifier_time: &mut Duration, msg: LayerReductionMessage<F>, layer: usize) -> (Vec<F>, F) {
        let inst = Instant::now();

        self.verifier.confirm_last_sum_check_msg(&msg, verifier, layer);

        let (next_layer_gate, claim) = verifier.handle_layer_reduction_message(&msg, s_i_plus_1);
        self.verifier.set_gate(next_layer_gate.as_slice(), layer + 1);
        *verifier_time += inst.elapsed();
        (next_layer_gate, claim)
    }

    fn track_creating_layer_reduc_msg<T: SumCheckProver<F>>(prover: &mut T, s_i_plus_1: usize, prover_time: &mut Duration) -> LayerReductionMessage<F> {
        let inst = Instant::now();
        let msg = prover.layer_reduction_message(s_i_plus_1);
        let elapsed = inst.elapsed();
        *prover_time += elapsed;
        println!("Time taken to create layer reduction msg: {:?}", elapsed);
        msg
    }

    pub fn run_circuit<T: SumCheckProver<F>> (
        &mut self,
    ) -> AnalysisResult
    {
        let mut an_res = AnalysisResult::new();
        let mut layer_connection = LayerConnection::new(vec![F::zero()], F::zero());
        // This is for all rounds except the last round.
        for i in 0..self.circuit.layers.len() {
            let (layer_conn, analysis) = self.handle_sum_check_for_layer::<T>(layer_connection, i);
            layer_connection = layer_conn;
            an_res.add_time_per_layer(analysis);
        }

        // This function panics == Verifier rejects.
        println!("Checking Final claim");

        let mut track = Track::new();
        let inst = Instant::now();
        self.verifier.verify_final_claimed_value_point(layer_connection.next_gate, layer_connection.claim_mi);
        let elapsed = inst.elapsed();
        track.add_verifier_time(elapsed);
        an_res.add_time_per_layer(track);

        an_res
    }

    fn handle_sum_check_for_layer<T: SumCheckProver<F>>(
        &mut self, mut layer_connection: LayerConnection<F>,
        i: usize
    ) -> (LayerConnection<F>, Track) {
        let s_i_plus_1 = self.get_correct_next_layer_size(i + 1);

        let value_extension = self.get_correct_value_extension(i + 1);
        if i == 0 {
            println!("Running initial round");
            self.handle_first_round::<T>(
                s_i_plus_1,
                &value_extension
            )
        } else {
            println!("Handling round {i}.");
            self.handle_intermediate_rounds::<T>(
                &mut layer_connection,
                s_i_plus_1,
                &value_extension,
                i,
            )
        }
    }

    fn get_correct_next_layer_size(&mut self, i: usize) -> usize {
        if i < self.circuit.layers.len() {
            let layer = &self.circuit.layers[i];
            let gates_len = layer.gates.len();
            log2_pow2(gates_len)
        } else {
            log2_pow2(self.input_layer.values.len())
        }
    }

    fn get_correct_value_extension(&mut self, i: usize) -> DenseMultilinearExtension<F> {
        if i < self.circuit.layers.len() {
            self.gkrprover
                .eval_circuit()
                .layers[i]
                .value_extension()
        } else {
            self.input_layer
                .value_extension()
        }
    }

    fn handle_intermediate_rounds<T: SumCheckProver<F>>(&mut self,
                                  layer_connection: &mut LayerConnection<F>,
                                  s_i_plus_1: usize,
                                  value_extension: &DenseMultilinearExtension<F>,
                                  layer: usize,
    ) -> (LayerConnection<F>, Track) {
        let (add_pred, mult_pred) = &self.gkrprover.predicates()[layer];
        let gkr_round: GKRRound<F> = GKRRound::new(&mult_pred.pred, &add_pred.pred, &value_extension, &value_extension);


        let (prover, elapsed_prover) = create_prover::<T, F>(&mut layer_connection.next_gate, gkr_round);

        let verifier =  StandardVerifier::new(2, layer_connection.claim_mi);
        let mut res = self.run_layer(prover, verifier, s_i_plus_1, layer);
        res.1.add_prover_time(elapsed_prover);
        res
    }

    fn handle_first_round<T: SumCheckProver<F>>(&mut self,
                          s_i_plus_1: usize,
                          value_extension: &DenseMultilinearExtension<F>,
    ) -> (LayerConnection<F>, Track)  {
        let output_claim = self.gkrprover.get_output_claim();
        let (add_pred, mult_pred) = &self.gkrprover.predicates()[0];
        let gkr_round: GKRRound<F> = GKRRound::new(&mult_pred.pred, &add_pred.pred, &value_extension.clone(), &value_extension.clone());

        // send claimed output to prover and get random gate for first iteration of sum-check.
        let (mut next_gate, mut elapsed_verifier) = self.get_random_gate_send_claim_to_verifier(&output_claim, 0);

        let (prover, elapsed_prover) = create_prover::<T, F>(&mut next_gate, gkr_round);

        let (verifier, time) = Self::create_initial_verifier_sumcheck(output_claim, &mut next_gate);
        elapsed_verifier += time;
        let mut res = self.run_layer(prover, verifier, s_i_plus_1, 0);
        res.1.add_prover_time(elapsed_prover);
        res.1.add_verifier_time(elapsed_verifier);
        res
    }

    /// This is only for the initial sum-check because it relies and sets the claim of the verifier
    /// based on the output claim sent by the prover.
    fn create_initial_verifier_sumcheck(output_claim: SparseMultilinearExtension<F>, next_gate: &mut Vec<F>) -> (StandardVerifier<F>, Duration) {
        let time = Instant::now();
        let verifier = StandardVerifier::new(2, output_claim.evaluate(&next_gate));
        let time = time.elapsed();
        (verifier, time)
    }

    fn get_random_gate_send_claim_to_verifier(&mut self, output_claim: &SparseMultilinearExtension<F>, layer: usize) -> (Vec<F>, Duration) {
        let instant = Instant::now();
        let initial_layer_size = self.get_correct_next_layer_size(0);
        let next_gate = self.verifier.random_gate(&output_claim, initial_layer_size);
        self.verifier.set_gate(next_gate.as_slice(), layer);
        let elapsed_verifier = instant.elapsed();
        (next_gate, elapsed_verifier)
    }

}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_ff::Field;
    use ark_poly::Polynomial;
    use ark_std::test_rng;
    use crate::gkr::gkr_circuit::compute_predicates;
    use crate::gkr::gkr_driver::GKRDriver;
    use crate::gkr::gkr_prover::GKRProver;
    use crate::gkr::gkr_round::GKRRound;
    use crate::gkr::gkr_verifier::GKRVerifier;
    use crate::gkr::layer::InputLayer;
    use crate::provers::fast::FastProver;
    use crate::structures::circuit_structures::GKRCircuit;
    use crate::structures::data_structures::{SumCheckProver, SumCheckVerifier};
    use crate::verifiers::standard_verifier::StandardVerifier;

    impl<F: Field> GKRDriver<F> {
        pub(crate) fn set_verifier(&mut self, p0: GKRVerifier<F>) {
            self.verifier = p0;
        }
    }

    #[test]
    fn simulate_full_gkr_round_with_naive_prover() {
        let layers = &[2, 4, 8, 32, 64, 128, 256, 512, 2048, 1024];
        let random_circuit: GKRCircuit<Fr> = GKRCircuit::random(layers, &mut test_rng());
        let input_layer: InputLayer<Fr> = InputLayer::random(layers.last().unwrap());

        let mut gkr_prover = GKRProver::new(random_circuit.clone(), input_layer.clone());

        let pred = compute_predicates(gkr_prover.circuit(), gkr_prover.input());
        gkr_prover.set_predicates(pred);
        let mut gkr_verifier = GKRVerifier::new(random_circuit.clone(), input_layer.clone());
        let pred = compute_predicates(gkr_prover.circuit(), gkr_prover.input());
        gkr_verifier.set_predicate(pred);
        let mut gkrdriver: GKRDriver<Fr> = GKRDriver::new(&gkr_prover, &gkr_verifier, random_circuit, input_layer);

        let gate = gkr_verifier.random_gate(&gkr_prover.get_output_claim(), 1);
        gkr_verifier.set_gate(&*gate, 0);
        gkrdriver.set_verifier(gkr_verifier);
        let pred = &compute_predicates(gkr_prover.circuit(), gkr_prover.input())[0];
        let vle = gkr_prover.eval_circuit().layers[1].value_extension();
        let gkr_round = GKRRound::new(&pred.1.pred, &pred.0.pred, &vle, &vle);
        let mut prover = FastProver::new(gkr_round.clone(), &gate);
        let verifier = StandardVerifier::new(2, prover.compute_sum());

        let (layer, _) = gkrdriver.run_layer(prover, verifier, vle.num_vars, 0);

        // Folded claim must match direct W_{i+1}(r_next) evaluation.
        let expected = gkr_round.vi().evaluate(&layer.next_gate);
        assert_eq!(layer.claim_mi, expected);
    }


}