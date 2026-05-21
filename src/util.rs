use std::iter;
use ark_ff::Field;
use ark_poly::{univariate, DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};
use ark_std::test_rng;

/// Originally taken from Arkworks.
/// Can be seen [here](https://github.com/arkworks-rs/sumcheck/blob/master/src/gkr_round_sumcheck/test.rs)
/// This should only really be used for testing. 
/// This has been modified to return two predicate functions for both addition and multiplication.
/// This is in order to support circuits that relies on both addition and multiplication gates.
pub fn random_gkr_round_gates<F: Field>(
    dim: usize,
) -> (
    SparseMultilinearExtension<F>,
    SparseMultilinearExtension<F>,
    DenseMultilinearExtension<F>,
    DenseMultilinearExtension<F>,
) {
    let rng = &mut test_rng();
    (
        SparseMultilinearExtension::rand_with_config(dim * 3, 1 << dim, rng),
        SparseMultilinearExtension::rand_with_config(dim * 3, 1 << dim, rng),
        DenseMultilinearExtension::rand(dim, rng),
        DenseMultilinearExtension::rand(dim, rng),
    )
}

/// Also taken from [arkworks](https://github.com/arkworks-rs/sumcheck/blob/master/src/gkr_round_sumcheck/test.rs)
/// Takes as input ${0,1}^n$ and returns $\mathbb{F}^n$
pub fn index_to_field_element<F: Field>(mut index: usize, mut nv: usize) -> Vec<F> {
    let mut ans = Vec::with_capacity(nv);
    while nv != 0 {
        ans.push(((index & 1) as u64).into());
        index >>= 1;
        nv -= 1;
    }
    ans
}

/// Again this is a testing function since it relies on ark_std::test_rng()
pub fn random_gate<F: Field>(label_length: usize) -> Vec<F> {
    let mut rng = test_rng();
    let mut res = Vec::with_capacity(label_length);
    (0..label_length).for_each(|_| {
        res.push(F::rand(&mut rng));
    });
    res
}


/// Takes two vectors representing the values picked for b*, c* in thalers 2025 opionated paper.
///
/// # Arguments
///
/// * `b_star`: The first fixed gate.
/// * `c_star`: The second fixed gate.
/// * `t`: The point on the line.
///
/// returns: Vec<F, Global>
/// Idk global has to be there otherwise it doesn't show F.
pub fn line_point<F: Field>(b_star: &[F], c_star: &[F], t: F) -> Vec<F> {
    assert_eq!(b_star.len(), c_star.len());
    b_star
        .iter()
        .zip(c_star.iter())
        .map(|(b, c)| {
            *b + t * (*c - *b)
        })
        .collect()
}


/// This function is taken from https://montekki.github.io/thaler-ch4-4/
pub fn line<F: Field>(b: &[F], c: &[F]) -> Vec<univariate::SparsePolynomial<F>> {
    iter::zip(b, c)
        .map(|(b, c)| {
            univariate::SparsePolynomial::from_coefficients_slice(&[(0, *b), (1, *c - b)])
        })
        .collect()
}