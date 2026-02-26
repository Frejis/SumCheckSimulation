use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};
use ark_std::test_rng;

/// Taken from arkworks sumcheck protocol.
/// Can be seen [here](https://github.com/arkworks-rs/sumcheck/blob/master/src/gkr_round_sumcheck/test.rs)
/// This should only really be used for testing.
pub fn random_gkr_round_gates<F: Field>(
    dim: usize,
) -> (
    SparseMultilinearExtension<F>,
    DenseMultilinearExtension<F>,
    DenseMultilinearExtension<F>,
) {
    let rng = &mut test_rng();
    (
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


pub fn line_point<F: Field>(b_star: &[F], c_star: &[F], t: F) -> Vec<F> {
    b_star
        .iter()
        .zip(c_star.iter())
        .map(|(b, c)| *b + t * (*c - *b))
        .collect()
}

pub fn poly_add_assign<F: Field>(lhs: &mut Vec<F>, rhs: &[F]) {
    if lhs.len() < rhs.len() {
        lhs.resize(rhs.len(), F::zero());
    }
    for (i, coeff) in rhs.iter().enumerate() {
        lhs[i] += *coeff;
    }
}

/// Multiply polynomial by (c0 + c1*x), coefficient form.
pub fn poly_mul_linear<F: Field>(poly: &[F], c0: F, c1: F) -> Vec<F> {
    let mut out = vec![F::zero(); poly.len() + 1];
    for (i, a) in poly.iter().enumerate() {
        out[i] += *a * c0;
        out[i + 1] += *a * c1;
    }
    out
}

/// Interpolate polynomial coefficients from samples (x_i, y_i) using Lagrange basis.
pub fn lagrange_interpolate_coeffs<F: Field>(xs: &[F], ys: &[F]) -> Vec<F> {
    assert_eq!(xs.len(), ys.len());
    let n = xs.len();
    let mut result = vec![F::zero(); n];

    for i in 0..n {
        // basis numerator: prod_{j != i} (x - x_j)
        let mut basis = vec![F::one()];
        let mut denom = F::one();

        for j in 0..n {
            if i == j {
                continue;
            }
            basis = poly_mul_linear(&basis, -xs[j], F::one());
            denom *= xs[i] - xs[j];
        }

        let inv = denom.inverse().expect("distinct interpolation points required");
        let scale = ys[i] * inv;
        for coeff in basis.iter_mut() {
            *coeff *= scale;
        }
        poly_add_assign(&mut result, &basis);
    }

    result
}

pub fn eval_univariate<F: Field>(coeffs: &[F], x: F) -> F {
    coeffs.iter().rev().fold(F::zero(), |acc, c| acc * x + c)
}