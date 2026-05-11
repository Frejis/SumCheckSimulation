use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, DenseUVPolynomial, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use ark_poly::univariate::DensePolynomial;
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
    for _ in 0..label_length {
        res.push(F::rand(&mut rng));
    }
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
pub fn _line_point<F: Field>(b_star: &[F], c_star: &[F], t: F) -> Vec<F> {
    assert_eq!(b_star.len(), c_star.len());
    b_star
        .iter()
        .zip(c_star.iter())
        .map(|(b, c)| {
            *b + t * (*c - *b)
        })
        .collect()
}


pub fn restrict_mle_to_line<F: Field>(
    mle: &DenseMultilinearExtension<F>,
    a: &[F],
    b: &[F],
    ts: &[F],
) -> Vec<F> {
    assert_eq!(a.len(), b.len());
    assert_eq!(a.len(), mle.num_vars);

    ts.iter()
        .map(|t| {
            let point: Vec<F> = a.iter()
                .zip(b.iter())
                .map(|(ai, bi)| *ai + (*bi - *ai) * *t)
                .collect();

            mle.evaluate(&point)
        })
        .collect()
}

pub fn interpolate_univariate<F: Field>(evals: &[F], points: &[F]) -> DensePolynomial<F> {
    assert_eq!(evals.len(), points.len());
    let mut poly = DensePolynomial::from_coefficients_vec(vec![F::zero()]);
    for i in 0..evals.len() {
        let mut term = DensePolynomial::from_coefficients_vec(vec![evals[i]]);
        for j in 0..evals.len() {
            if i == j {
                continue;
            }
            let denominator = (points[i] - points[j]).inverse().unwrap();
            // (X - points[j]) * denominator = denominator * X - points[j] * denominator
            let sub_poly =
                DensePolynomial::from_coefficients_vec(vec![-points[j] * denominator, denominator]);
            term = term.naive_mul(&sub_poly);
        }
        poly = poly + term;
    }
    poly
}

#[cfg(test)]
mod test {
    use ark_bls12_381::Fr;
    use ark_ff::{One, Zero};
    use crate::util::{_line_point};

    #[test]
    fn test_line_point() {
        let size = 5;
        let b_star = vec![Fr::zero(); size];
        let c_star = vec![Fr::one(); size];
        assert_eq!(_line_point(b_star.as_ref(), c_star.as_ref(), Fr::zero()), b_star);
        assert_eq!(_line_point(b_star.as_ref(), c_star.as_ref(), Fr::one()), c_star)
    }

}