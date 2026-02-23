use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};
use ark_std::rand::RngCore;
use ark_std::test_rng;

/// Taken from arkworks sumcheck protocol.
/// Can be seen [here](https://github.com/arkworks-rs/sumcheck/blob/master/src/gkr_round_sumcheck/test.rs)
/// This should only really be used for testing.
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
    /*
    let mut res = Vec::with_capacity(label_length);
    for i in 0..label_length {
        res[i] = F::rand(&mut rng);
    }
    res
    */
    let mut res = Vec::with_capacity(label_length);
    for _ in 0..label_length {
        res.push(F::rand(&mut rng));
    }
    res
}