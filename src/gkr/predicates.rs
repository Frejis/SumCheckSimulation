use ark_ff::Field;
use ark_poly::SparseMultilinearExtension;

pub struct MultPredicate<F: Field> {
    pub pred: SparseMultilinearExtension<F>,
}

impl<F: Field> MultPredicate<F> {
    pub fn new(pred: SparseMultilinearExtension<F>) -> Self {
        Self { pred }
    }
}

pub struct AddPredicate<F: Field> {
    pub pred: SparseMultilinearExtension<F>,
}

impl<F: Field> AddPredicate<F> {
    pub fn new(pred: SparseMultilinearExtension<F>) -> Self {
        Self { pred }
    }
}