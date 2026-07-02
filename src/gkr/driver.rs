//! The GKR driver: simulates the full protocol between a [`GKRProver`] and a
//! [`GKRVerifier`], running one sum-check per circuit layer and timing prover
//! and verifier work separately.

use ark_ff::Field;
use ark_poly::univariate::SparsePolynomial;
use ark_poly::{DenseMultilinearExtension, Polynomial};

use crate::gkr::{GKRCircuit, GKRProver, GKRVerifier, InputLayer};
use crate::sumcheck::{GKRRound, StandardVerifier, SumCheckProver, SumCheckVerifier};
use crate::timing::{timed, AnalysisResult, Track};
use crate::util::log2_pow2;

/// Degree bound handed to the sum-check verifier. The round polynomials sent
/// in this simulation are degree 1 (see
/// [`SumCheckProver::get_verifier_function`]).
const MAX_ROUND_POLY_DEGREE: usize = 2;

/// Connects two consecutive layers of the protocol: the random point the
/// verifier picked during layer `i`'s reduction and the claimed value of
/// `W_{i+1}` at that point.
pub struct LayerConnection<F: Field> {
    pub next_gate: Vec<F>,
    pub claim_mi: F,
}

impl<F: Field> LayerConnection<F> {
    pub fn new(next_gate: Vec<F>, claim_mi: F) -> Self {
        Self { next_gate, claim_mi }
    }
}

/// Runs the GKR protocol over a circuit, generic over the sum-check prover
/// used for each layer.
pub struct GKRDriver<F: Field> {
    prover: GKRProver<F>,
    verifier: GKRVerifier<F>,
    circuit: GKRCircuit<F>,
    input_layer: InputLayer<F>,
}

impl<F: Field> GKRDriver<F> {
    pub fn new(
        prover: GKRProver<F>,
        verifier: GKRVerifier<F>,
        circuit: GKRCircuit<F>,
        input_layer: InputLayer<F>,
    ) -> Self {
        Self { prover, verifier, circuit, input_layer }
    }

    /// Simulates the whole protocol: one sum-check per layer followed by the
    /// final input-layer check. Panicking == the verifier rejects.
    pub fn run_circuit<P: SumCheckProver<F>>(&mut self) -> AnalysisResult {
        let mut results = AnalysisResult::new();
        let mut connection = LayerConnection::new(vec![F::zero()], F::zero());

        for layer in 0..self.circuit.layers.len() {
            let (next_connection, track) = self.run_layer_sum_check::<P>(connection, layer);
            connection = next_connection;
            results.add_time_per_layer(track);
        }

        let mut track = Track::new();
        let (_, elapsed) = timed(|| {
            self.verifier
                .verify_final_claimed_value_point(&connection.next_gate, connection.claim_mi)
        });
        track.add_verifier_time(elapsed);
        results.add_time_per_layer(track);

        results
    }

    /// Runs the sum-check for one layer (2·s_{i+1} rounds) plus the layer
    /// reduction, and returns the connection to the next layer.
    fn run_layer_sum_check<P: SumCheckProver<F>>(
        &mut self,
        connection: LayerConnection<F>,
        layer: usize,
    ) -> (LayerConnection<F>, Track) {
        let s_i_plus_1 = self.layer_num_vars(layer + 1);
        let value_extension = self.layer_value_extension(layer + 1);
        let round = {
            let predicates = &self.prover.predicates()[layer];
            GKRRound::new(&predicates.mult, &predicates.add, &value_extension, &value_extension)
        };

        let mut track = Track::new();

        let (gate, claim) = if layer == 0 {
            // The prover opens with the claimed output MLE; the verifier picks
            // a random output gate g and the sum-check starts on ~W_0(g).
            let output_claim = self.prover.output_claim();
            let output_vars = self.layer_num_vars(0);
            let ((gate, claim), elapsed) = timed(|| {
                let gate = self.verifier.random_gate(output_vars);
                self.verifier.set_gate(&gate, 0);
                let claim = output_claim.evaluate(&gate);
                (gate, claim)
            });
            track.add_verifier_time(elapsed);
            (gate, claim)
        } else {
            (connection.next_gate, connection.claim_mi)
        };

        let (mut prover, elapsed) = timed(|| P::new(round, &gate));
        track.add_prover_time(elapsed);
        let mut verifier = StandardVerifier::new(MAX_ROUND_POLY_DEGREE, claim);

        for _ in 0..2 * s_i_plus_1 {
            let (points, elapsed) = timed(|| prover.get_verifier_function());
            track.add_prover_time(elapsed);

            // r_j is the verifier's random point for this round.
            let (r_j, elapsed) =
                timed(|| verifier.handle_round(&SparsePolynomial::from_coefficients_vec(points)));
            track.add_verifier_time(elapsed);

            let (new_claim, elapsed) = timed(|| {
                prover.fix_variable(r_j);
                prover.compute_sum()
            });
            track.add_prover_time(elapsed);

            let (_, elapsed) = timed(|| verifier.set_claim(new_claim));
            track.add_verifier_time(elapsed);
        }

        let (msg, elapsed) = timed(|| prover.layer_reduction_message(s_i_plus_1));
        track.add_prover_time(elapsed);

        let ((next_gate, next_claim), elapsed) = timed(|| {
            self.verifier.confirm_last_sum_check_msg(&msg, &verifier, layer);
            let (next_gate, next_claim) = verifier.handle_layer_reduction_message(&msg, s_i_plus_1);
            self.verifier.set_gate(&next_gate, layer + 1);
            (next_gate, next_claim)
        });
        track.add_verifier_time(elapsed);

        (LayerConnection::new(next_gate, next_claim), track)
    }

    /// Number of gate-label variables of `layer`; past the last circuit layer
    /// this is the input layer's.
    fn layer_num_vars(&self, layer: usize) -> usize {
        if layer < self.circuit.layers.len() {
            log2_pow2(self.circuit.layers[layer].gates.len())
        } else {
            log2_pow2(self.input_layer.values.len())
        }
    }

    /// The value MLE `~W_layer`; past the last circuit layer this is the input
    /// layer's extension.
    fn layer_value_extension(&self, layer: usize) -> DenseMultilinearExtension<F> {
        if layer < self.circuit.layers.len() {
            self.prover.evaluated_circuit().layers[layer].value_extension()
        } else {
            self.input_layer.value_extension()
        }
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_ff::Zero;
    use ark_poly::Polynomial;
    use ark_std::test_rng;

    use crate::gkr::{compute_predicates, GKRCircuit, GKRProver, GKRVerifier, InputLayer};
    use crate::sumcheck::{FastProver, NaiveProver, SumCheckProver};

    use super::{GKRDriver, LayerConnection};

    fn driver(layer_sizes: &[usize]) -> GKRDriver<Fr> {
        let circuit: GKRCircuit<Fr> = GKRCircuit::random(layer_sizes, &mut test_rng());
        let input: InputLayer<Fr> = InputLayer::random(*layer_sizes.last().unwrap());
        let predicates = compute_predicates(&circuit, &input);
        let prover = GKRProver::new(&circuit, &input, predicates.clone());
        let verifier = GKRVerifier::new(input.clone(), predicates);
        GKRDriver::new(prover, verifier, circuit, input)
    }

    #[test]
    fn layer_sum_check_folds_claim_to_next_layer_evaluation() {
        let mut driver = driver(&[2, 4, 8, 32, 64, 128, 256, 512, 2048, 1024]);

        let (connection, _) = driver
            .run_layer_sum_check::<FastProver<Fr>>(LayerConnection::new(Vec::new(), Fr::zero()), 0);

        // The folded claim must match a direct evaluation of W_1 at the point
        // chosen by the verifier.
        let w_1 = driver.prover.evaluated_circuit().layers[1].value_extension();
        assert_eq!(connection.claim_mi, w_1.evaluate(&connection.next_gate));
    }

    fn assert_full_protocol_accepts_honest_prover<P: SumCheckProver<Fr>>() {
        // run_circuit panics if the verifier rejects at any point.
        driver(&[2, 4, 8, 16]).run_circuit::<P>();
    }

    #[test]
    fn full_protocol_accepts_honest_naive_prover() {
        assert_full_protocol_accepts_honest_prover::<NaiveProver<Fr>>();
    }

    #[test]
    fn full_protocol_accepts_honest_fast_prover() {
        assert_full_protocol_accepts_honest_prover::<FastProver<Fr>>();
    }
}
