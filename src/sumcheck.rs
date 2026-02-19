// These structs are taken from ARKWORKS
// They are intended to help defining a GKRROUND.
// THe purpose is for this to help align my interface with the implementation of ARKWORKS
// This is done to enhance benchmarking between different implementation of sum-check.


use std::marker::PhantomData;
use ark_ff::Field;
use ark_poly::{
    DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension
};
use ark_std::rand::RngCore;
use ark_test_curves::ark_ff;
use crate::data_structures::GKRRoundSumcheckSubClaim;
use crate::rng::FeedableRNG;
use crate::errors::Error;

// Taken from ARKWORKs modfied derive tratis.
/// Prover Message
#[derive(Clone)]
pub struct ProverMsg<F: Field> {
    /// evaluations on P(0), P(1), P(2), ...
    pub(crate) evaluations: Vec<F>,
}

/// Sumcheck Argument for GKR Round function
pub struct GKRRoundSumcheck<F: Field> {
    _marker: PhantomData<F>,
}

// This GKRProof is taken from the ARKWORKS definition.
pub struct GKRProof<F: Field> {
    pub(crate) sumcheck_msgs: Vec<ProverMsg<F>>,
}

impl<F: Field> GKRProof<F> {
    /// Extract the witness (i.e. the sum of GKR)
    pub fn extract_sum(&self) -> F {
        todo!()
    }
}

pub fn random_gkr_instance<F: Field, R: RngCore>(
    dim: usize,
    rng: &mut R,
) -> (
    SparseMultilinearExtension<F>,
    DenseMultilinearExtension<F>,
    DenseMultilinearExtension<F>,
) {
    (
        SparseMultilinearExtension::rand_with_config(dim * 3, 1 << dim, rng),
        DenseMultilinearExtension::rand(dim, rng),
        DenseMultilinearExtension::rand(dim, rng),
    )
}

pub trait SumcheckProver<F: ark_ff::Field> {
    fn compute_sum(
        f1: &SparseMultilinearExtension<F>,
        f2: &DenseMultilinearExtension<F>,
        f3: &DenseMultilinearExtension<F>,
        r: &[F],
    ) -> F;
}

// The interface is taken from ARKWORKS, however i have added comments and "rewritten it" such that it is my solution.
impl<F: Field> GKRRoundSumcheck<F> {
    /// Takes a GKR Round Function and input, prove the sum.
    /// * `f1`,`f2`,`f3`: represents the GKR round function
    /// * `g`: represents the fixed input.
    pub fn prove<R: FeedableRNG>(
        rng: &mut R,
        f1: &SparseMultilinearExtension<F>,
        f2: &DenseMultilinearExtension<F>,
        f3: &DenseMultilinearExtension<F>,
        g: &[F],
    ) -> GKRProof<F> {
        // In GKR protocol we have mult_k(r, i, j).
        // Hence the amount of variables will be 3 * variables in V_(k-1)(i) and V_(k-1)(j).
        assert_eq!(f1.num_vars, 3 * f2.num_vars);
        assert_eq!(f1.num_vars, 3 * f3.num_vars);
    
        
        todo!()
    }

    /// Takes a GKR Round Function, input, and proof, and returns a subclaim.
    ///
    /// If the `claimed_sum` is correct, then it is `subclaim.verify_subclaim` will return true.
    /// Otherwise, it is very likely that `subclaim.verify_subclaim` will return false.
    /// Larger field size guarantees smaller soundness error.
    /// * `f2_num_vars`: represents number of variables of f2
    pub fn verify<R: FeedableRNG>(
        rng: &mut R,
        f2_num_vars: usize,
        proof: &GKRProof<F>,
        claimed_sum: F,
    ) -> Result<GKRRoundSumcheckSubClaim<F>, Error> {
        todo!()
    }
}