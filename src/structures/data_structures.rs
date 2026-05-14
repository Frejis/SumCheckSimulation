/*
This file will contain structures relevant to setting up the proof system.
*/
use std::iter;
use ark_ff::{Field, Zero};
use ark_poly::{univariate, DenseUVPolynomial, MultilinearExtension, SparseMultilinearExtension};
use crate::gkr::gkr_round::GKRRound;
use crate::gkr::layer::LayerReductionMessage;

pub trait SumCheckProver<F: Field> {
    // Computes the sum so we can have an alleged claim of the functions.
    fn compute_sum(&mut self) -> F;

    // Creates a function that has one variable (meaning it fixes all other variables)
    fn get_verifier_function(&mut self) -> SparseMultilinearExtension<F>;

    fn fix_variable(&mut self, random_field_element: F);

    fn layer_reduction_message(&self, s_i_plus_1: usize) -> LayerReductionMessage<F>;

    /// This function is taken from: https://montekki.github.io/thaler-ch4-4/
    fn restrict_poly<M: MultilinearExtension<F>>(
        b: &[F],
        c: &[F],
        mle: &M,
    ) -> univariate::SparsePolynomial<F> {
        let k: Vec<_> = iter::zip(b, c).map(|(b, c)| *c - b).collect();

        let evaluations = mle.to_evaluations();
        let num_vars = mle.num_vars();

        let mut res = univariate::SparsePolynomial::zero();
        

        for (i, evaluation) in evaluations.iter().enumerate() {
            let mut p = univariate::SparsePolynomial::from_coefficients_vec(vec![(0, *evaluation)]);
            for bit in 0..num_vars {
                let mut b =
                    univariate::SparsePolynomial::from_coefficients_vec(vec![(0, b[bit]), (1, k[bit])]);

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
    fn verify_degree(&self, fx: &SparseMultilinearExtension<F>) -> bool;

    // Returns a random field element from the verifier
    fn get_random_field_element(&mut self) -> F;

    // Takes as input a multilinear extension and checks that for each field their sum is the claim.
    fn check_claimed_value(&self, fx: &SparseMultilinearExtension<F>) -> bool;

    /// Should ideally take a function by the prover and do all necessary checks
    /// If any fails then it panics, and if everything is good then it returns a random field element.
    fn handle_round(&mut self, fx: &SparseMultilinearExtension<F>) -> F;

    fn set_claim(&mut self, claim: F);

    fn final_check(&self);
}