/*
This file will contain structures relevant to setting up the proof system.
*/
use std::iter;
use std::time::Duration;
use ark_ff::{Field, Zero};
use ark_poly::{univariate, DenseMultilinearExtension, DenseUVPolynomial, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use ark_poly::univariate::SparsePolynomial;
use crate::gkr::gkr_round::GKRRound;
use crate::gkr::layer::LayerReductionMessage;
use crate::gkr::predicates::{AddPredicate, MultPredicate};

pub trait SumCheckProver<F: Field> {
    // Computes the sum so we can have an alleged claim of the functions.
    fn compute_sum(&mut self) -> F;

    // Creates a function that has one variable (meaning it fixes all other variables)
    fn get_verifier_function(&mut self) -> Vec<(usize, F)>;

    fn fix_variable(&mut self, random_field_element: F);

    fn layer_reduction_message(&self, s_i_plus_1: usize) -> LayerReductionMessage<F>;

    /// This function is taken from: https://montekki.github.io/thaler-ch4-4/
    fn restrict_poly<M: MultilinearExtension<F>>(
        b: &[F],
        c: &[F],
        mle: &M,
    ) -> SparsePolynomial<F> {
        let k: Vec<_> = iter::zip(b, c).map(|(b, c)| *c - b).collect();

        let evaluations = mle.to_evaluations();
        let num_vars = mle.num_vars();
        
        let mut res = SparsePolynomial::zero();


        for (i, evaluation) in evaluations.iter().enumerate() {
            let mut p = SparsePolynomial::from_coefficients_vec(vec![(0, *evaluation)]);
            for bit in 0..num_vars {
                let mut b =
                    SparsePolynomial::from_coefficients_vec(vec![(0, b[bit]), (1, k[bit])]);

                if i & (1 << bit) == 0 {
                    b = (&univariate::DensePolynomial::from_coefficients_vec(vec![F::one()]) - &b)
                        .into();
                }

                p = p.mul(&b);
            }

            res += &p;
        }

        res
    }

    fn new(gkr_round: GKRRound<F>, gate: &[F]) -> Self;
}

pub trait SumCheckVerifier<F: Field> {
    // Has to check the degree of the function to ensure no one cheats.
    fn verify_degree(&self, fx: &SparsePolynomial<F>) -> bool;

    // Returns a random field element from the verifier
    fn get_random_field_element(&mut self) -> F;

    // Takes as input a multilinear extension and checks that for each field their sum is the claim.
    fn check_claimed_value(&self, fx: &SparsePolynomial<F>) -> bool;

    /// Should ideally take a function by the prover and do all necessary checks
    /// If any fails then it panics, and if everything is good then it returns a random field element.
    fn handle_round(&mut self, fx: &SparsePolynomial<F>) -> F;

    fn set_claim(&mut self, claim: F);


    /// This should only be used during Sum-check simulation as the final check is otherwise
    /// handled by checking that the layer is correct. As such this is given the MLE as input.
    fn final_check(&self, gate: &[F], add_pred: &SparseMultilinearExtension<F>, mult_pred: &SparseMultilinearExtension<F>, mle: DenseMultilinearExtension<F>, vrf_func: SparsePolynomial<F>);
    
    fn new(max_degree: usize, claim: F) -> Self;
}

pub struct AnalysisResult {
    time_per_layer: Vec<Track>
}

pub struct Track {
    prover: Duration,
    verifier: Duration,
}

impl Track {
    pub fn new() -> Self {
        Self {
            prover: Duration::new(0,0),
            verifier: Duration::new(0,0),
        }
    }

    pub fn new_times(prover: Duration, verifier: Duration) -> Self {
        Self {
            prover,
            verifier,
        }
    }

    pub fn add_verifier_time(&mut self, time: Duration) {
        self.verifier += time;
    }

    pub fn add_prover_time(&mut self, time: Duration) {
        self.prover += time;
    }

    pub fn prover(&self) -> Duration {
        self.prover
    }

    pub fn verifier(&self) -> Duration {
        self.verifier
    }
}

impl AnalysisResult {
    pub fn new() -> Self {
        Self { time_per_layer: Vec::new() }
    }

    pub fn add_verifier_time(&mut self, time: Duration, layer: usize) {
        self.time_per_layer[layer].add_verifier_time(time);
    }

    pub fn add_time_per_layer(&mut self, prv_time: Track) {
        self.time_per_layer.push(prv_time)
    }

    pub fn print_total(&self) {
        let mut prover = Duration::new(0,0);
        let mut verifier = Duration::new(0,0);
        for i in self.time_per_layer.iter() {
            prover += i.prover;
            verifier += i.verifier;
        }
        println!("Prover spent: {:?}. Verifier spent: {:?}", prover, verifier);
    }

    pub fn print_each_layer(&self) {
        for (layer, track) in self.time_per_layer.iter().enumerate() {
            println!("Layer {layer}. Prover spent: {:?}. Verifier spent: {:?}.",
                track.prover,
                track.verifier
            )
        }
    }

    pub fn get_time_for_layer(&self, layer: usize) -> &Track {
        &self.time_per_layer[layer]
    }

    pub fn add_prover_time(&mut self, time: Duration, layer: usize) {
        self.time_per_layer[layer].add_prover_time(time);
    }
}

