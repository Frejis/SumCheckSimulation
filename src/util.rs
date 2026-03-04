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

/// Univariate Lagrange interpolation at point `r` using evaluations at 0, 1, ..., n-1.
pub fn interpolate_univariate<F: Field>(evals: &[F], r: F) -> F {
    let n = evals.len();
    let mut result = F::zero();
    let mut x_coords = Vec::with_capacity(n);
    for i in 0..n {
        x_coords.push(F::from(i as u64));
    }
    
    for j in 0..n {
        let mut num = F::one();
        let mut den = F::one();
        for i in 0..n {
            if i == j {
                continue;
            }
            num *= r - x_coords[i];
            den *= x_coords[j] - x_coords[i];
        }
        result += evals[j] * num * den.inverse().unwrap();
    }
    result
}

/// This is an implementation of lagrange interpolation on multilinear polynomials.
/// https://www.geeksforgeeks.org/dsa/lagranges-interpolation/
/// # Arguments
///
/// * `b`: &[F]
/// * `c`: &[F]
///
/// returns: Vec<F, Global>
/// The math in the definition is labeled to align with terminology in ProofsArgsAndZk.
/// Computes the Lagrange basis polynomials as defined in Eq 3.2.
/// 'a' is the evaluation point (x_1, ..., x_v).
/// Returns a vector of size 2^v where each element is chi_w(x).
pub fn _lagrange_interpolate_coeffs<F: Field>(a: &[F]) -> Vec<F> {
    let v = a.len();
    let mut basis = Vec::with_capacity(1 << v);

    // Start with the base case (empty product is 1)
    basis.push(F::one());

    // Iteratively apply Eq 3.2: chi_w(x) = product (x_i*w_i + (1-x_i)(1-w_i))
    for &xi in a {
        let mut next_basis = Vec::with_capacity(basis.len() * 2);
        for prev in basis {
            // Case w_i = 0: term is (1 - x_i)
            next_basis.push(prev * (F::one() - xi));
            // Case w_i = 1: term is x_i
            next_basis.push(prev * xi);
        }
        basis = next_basis;
    }

    basis
}

#[cfg(test)]
mod test {
    use ark_bls12_381::Fr;
    use ark_ff::{Field, One, Zero};
    use crate::util::{_line_point, _lagrange_interpolate_coeffs};

    #[test]
    fn test_line_point() {
        let size = 5;
        let b_star = vec![Fr::zero(); size];
        let c_star = vec![Fr::one(); size];
        assert_eq!(_line_point(b_star.as_ref(), c_star.as_ref(), Fr::zero()), b_star);
        assert_eq!(_line_point(b_star.as_ref(), c_star.as_ref(), Fr::one()), c_star)
    }

    #[test]
    fn test_lagrange_coeffs_boolean_property() {
        // Point (1, 0) corresponds to index 2 in lexicographical order (1*2^1 + 0*2^0)
        // because our implementation follows Equation 3.2's product structure.
        let point = vec![Fr::one(), Fr::zero()];
        let coeffs = _lagrange_interpolate_coeffs(&point);

        // For v=2, there are 2^2 = 4 coefficients
        assert_eq!(coeffs.len(), 4);

        // According to Lemma 3.6, at a boolean point w*, chi_w*(w*) = 1 and others are 0.
        // Index 2 is w=(1,0)
        assert_eq!(coeffs[0], Fr::zero()); // w=(0,0)
        assert_eq!(coeffs[1], Fr::zero()); // w=(0,1)
        assert_eq!(coeffs[2], Fr::one());  // w=(1,0) -> Correct!
        assert_eq!(coeffs[3], Fr::zero()); // w=(1,1)
    }

    #[test]
    fn test_lagrange_coeffs_summation() {
        // Test with arbitrary field elements (interpolation point outside the hypercube)
        let point = vec![Fr::from(5u64), Fr::from(12u64), Fr::from(7u64)];
        let coeffs = _lagrange_interpolate_coeffs(&point);

        // The sum of all chi_w(x) must always be 1
        let sum: Fr = coeffs.iter().sum();
        assert_eq!(sum, Fr::one());
    }

    #[test]
    fn test_specific_value() {
        // For v=1 at x=0.5, chi_0 = 0.5 and chi_1 = 0.5
        // Using Fr::from(2).inverse() to get 0.5 in the field
        let half = Fr::from(2u64).inverse().unwrap();
        let point = vec![half];
        let coeffs = _lagrange_interpolate_coeffs(&point);

        assert_eq!(coeffs[0], half); // (1 - 0.5)
        assert_eq!(coeffs[1], half); // (0.5)
    }
}