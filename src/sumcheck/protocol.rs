//! The prover/verifier interfaces of the sum-check protocol and the
//! layer-reduction message that connects consecutive GKR layers.

use std::iter;

use ark_ff::{Field, Zero};
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};
use ark_poly::{
    DenseMultilinearExtension, DenseUVPolynomial, MultilinearExtension, Polynomial,
    SparseMultilinearExtension,
};

use crate::sumcheck::GKRRound;

/// Message sent by the prover after the last sum-check round of a GKR layer.
///
/// It reduces the verifier's two claims about `W_{i+1}` — at the points `b*`
/// and `c*` fixed during sum-check — to a single claim about one point on the
/// line through `b*` and `c*`.
#[derive(Clone, Debug)]
pub struct LayerReductionMessage<F: Field> {
    z1: F,
    z2: F,
    qt: SparsePolynomial<F>,
}

impl<F: Field> LayerReductionMessage<F> {
    pub fn new(z1: F, z2: F, qt: SparsePolynomial<F>) -> Self {
        Self { z1, z2, qt }
    }

    /// The claimed value `W_{i+1}(b*)`.
    pub fn z1(&self) -> F {
        self.z1
    }

    /// The claimed value `W_{i+1}(c*)`.
    pub fn z2(&self) -> F {
        self.z2
    }

    /// The restriction `q(t) = W_{i+1}(l(t))` of the layer MLE to the line
    /// through `b*` and `c*`.
    pub fn qt(&self) -> &SparsePolynomial<F> {
        &self.qt
    }
}

/// A prover for one sum-check execution over a [`GKRRound`] instance.
pub trait SumCheckProver<F: Field> {
    /// Creates the prover for `round` with the gate label already fixed to `gate`.
    fn new(round: GKRRound<F>, gate: &[F]) -> Self;

    /// Computes the sum over the remaining boolean hypercube, i.e. the claim
    /// the prover currently stands behind.
    fn compute_sum(&mut self) -> F;

    /// Returns the current round polynomial as the pair
    /// `[(0, g(0)), (1, g(1))]`, where `g` fixes the next free variable.
    /// Summing the entries therefore gives `g(0) + g(1)`, which the verifier
    /// checks against the current claim.
    fn get_verifier_function(&mut self) -> Vec<(usize, F)>;

    /// Fixes the next free variable to the verifier's random point `r`.
    fn fix_variable(&mut self, r: F);

    /// All random points fixed so far, in order.
    fn fixed_variables(&self) -> &[F];

    /// The untouched MLE `W_{i+1}` of the next layer's values, used for the
    /// layer-reduction message.
    fn layer_value_mle(&self) -> &DenseMultilinearExtension<F>;

    /// Builds the layer-reduction message from the points fixed during this
    /// sum-check: `b*` are the first `s_i_plus_1` fixed variables, `c*` the
    /// next `s_i_plus_1`.
    fn layer_reduction_message(&self, s_i_plus_1: usize) -> LayerReductionMessage<F> {
        let mle = self.layer_value_mle();
        let b_star = self.fixed_variables()[..s_i_plus_1].to_vec();
        let c_star = self.fixed_variables()[s_i_plus_1..2 * s_i_plus_1].to_vec();
        let qt = restrict_poly(&b_star, &c_star, mle);
        let z1 = mle.evaluate(&b_star);
        let z2 = mle.evaluate(&c_star);
        LayerReductionMessage::new(z1, z2, qt)
    }
}

/// A verifier for one sum-check execution.
pub trait SumCheckVerifier<F: Field> {
    fn new(max_degree: usize, claim: F) -> Self;

    /// Checks the degree of the round polynomial to ensure the prover does not cheat.
    fn verify_degree(&self, fx: &SparsePolynomial<F>) -> bool;

    /// Samples (and records) a random field element to send to the prover.
    fn get_random_field_element(&mut self) -> F;

    /// Checks that the round polynomial is consistent with the current claim.
    fn check_claimed_value(&self, fx: &SparsePolynomial<F>) -> bool;

    /// Performs all checks for one round, panicking on failure, and returns
    /// the verifier's random point for this round.
    fn handle_round(&mut self, fx: &SparsePolynomial<F>) -> F;

    fn set_claim(&mut self, claim: F);

    /// Final oracle check of a standalone sum-check simulation. Inside GKR the
    /// final check is instead handled by reducing to a claim about the next
    /// layer, which is why this is given the layer MLE directly.
    fn final_check(
        &self,
        gate: &[F],
        add_pred: &SparseMultilinearExtension<F>,
        mult_pred: &SparseMultilinearExtension<F>,
        mle: DenseMultilinearExtension<F>,
        vrf_func: Vec<(usize, F)>,
    );
}

/// Restricts `mle` to the line through `b` and `c`, returning the univariate
/// polynomial `q(t) = mle(l(t))` with `l(0) = b` and `l(1) = c`.
///
/// Taken from <https://montekki.github.io/thaler-ch4-4/>.
pub fn restrict_poly<F: Field, M: MultilinearExtension<F>>(
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
            let mut term =
                SparsePolynomial::from_coefficients_vec(vec![(0, b[bit]), (1, k[bit])]);

            if i & (1 << bit) == 0 {
                term = (&DensePolynomial::from_coefficients_vec(vec![F::one()]) - &term).into();
            }

            p = p.mul(&term);
        }

        res += &p;
    }

    res
}
