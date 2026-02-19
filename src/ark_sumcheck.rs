use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};

fn index_to_field_element<F: Field>(mut index: usize, mut nv: usize) -> Vec<F> {
    let mut ans = Vec::with_capacity(nv);
    while nv != 0 {
        ans.push(((index & 1) as u64).into());
        index >>= 1;
        nv -= 1;
    }
    ans
}

pub fn calculate_sum_naive<F: Field>(
    f1: &SparseMultilinearExtension<F>,
    f2: &DenseMultilinearExtension<F>,
    f3: &DenseMultilinearExtension<F>,
    g: &[F],
) -> F {
    let dim = f2.num_vars;
    assert_eq!(f1.num_vars, 3 * dim);
    assert_eq!(f3.num_vars, dim);
    let f1_g = f1.fix_variables(g);
    let mut sum_xy = F::zero();
    for x in 0..(1 << dim) {
        let f2_x = f2[x];
        let f1_gx = f1_g
            .fix_variables(&index_to_field_element(x, dim))
            .to_dense_multilinear_extension();
        for y in 0..(1 << dim) {
            sum_xy += f1_gx[y] * f2_x * f3[y];
        }
    }
    sum_xy
}