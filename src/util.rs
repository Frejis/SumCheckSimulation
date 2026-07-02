//! Small helpers shared across the crate.

use ark_ff::Field;
use ark_std::test_rng;

/// Returns `log2(n)` for a power of two and panics otherwise.
pub fn log2_pow2(n: usize) -> usize {
    assert!(n.is_power_of_two());
    n.trailing_zeros() as usize
}

/// Maps an index in `{0,1}^nv` to the corresponding vector of field elements
/// (least significant bit first).
///
/// Taken from [arkworks](https://github.com/arkworks-rs/sumcheck/blob/master/src/gkr_round_sumcheck/test.rs).
pub fn index_to_field_element<F: Field>(mut index: usize, mut nv: usize) -> Vec<F> {
    let mut ans = Vec::with_capacity(nv);
    while nv != 0 {
        ans.push(((index & 1) as u64).into());
        index >>= 1;
        nv -= 1;
    }
    ans
}

/// Evaluates the unique line `l` with `l(0) = b_star` and `l(1) = c_star` at
/// the point `t` (see Thaler's 2025 opinionated survey).
pub fn line_point<F: Field>(b_star: &[F], c_star: &[F], t: F) -> Vec<F> {
    assert_eq!(b_star.len(), c_star.len());
    b_star
        .iter()
        .zip(c_star.iter())
        .map(|(b, c)| *b + t * (*c - *b))
        .collect()
}

/// Returns a random gate label of the given length.
///
/// Only for testing/simulation: relies on the deterministic `ark_std::test_rng`.
pub fn random_gate<F: Field>(label_length: usize) -> Vec<F> {
    let mut rng = test_rng();
    (0..label_length).map(|_| F::rand(&mut rng)).collect()
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fr;
    use ark_ff::{One, Zero};

    use super::{index_to_field_element, line_point};

    #[test]
    fn index_to_field_element_matches_bit_decomposition() {
        for mask in 0..30 {
            let field_index: Vec<Fr> = index_to_field_element(mask, 30);
            let expected: Vec<Fr> = (0..30)
                .map(|j| if (mask >> j) & 1 != 0 { Fr::one() } else { Fr::zero() })
                .collect();
            assert_eq!(field_index, expected);
        }
    }

    #[test]
    fn line_point_hits_endpoints() {
        let size = 5;
        let b_star = vec![Fr::zero(); size];
        let c_star = vec![Fr::one(); size];
        assert_eq!(line_point(&b_star, &c_star, Fr::zero()), b_star);
        assert_eq!(line_point(&b_star, &c_star, Fr::one()), c_star);
    }
}
