use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};

pub struct FastProver<F: Field> {
    mult: SparseMultilinearExtension<F>,
    fixed_mult: SparseMultilinearExtension<F>,
    gate: Vec<F>,
    vi: DenseMultilinearExtension<F>,
    vj: DenseMultilinearExtension<F>,
    p: Vec<F>,
    q: Vec<F>,
}

impl<F: Field> FastProver<F> {
    pub fn new(
        mult: SparseMultilinearExtension<F>,
        gate: Vec<F>,
        vi: DenseMultilinearExtension<F>,
        vj: DenseMultilinearExtension<F>,
    ) -> Self {
        Self {
            fixed_mult: mult.fix_variables(&*gate),
            mult,
            gate,
            vi,
            vj,
            p: vec![],
            q: vec![],
        }
    }


}