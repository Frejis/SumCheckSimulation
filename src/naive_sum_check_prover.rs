use std::ops::AddAssign;

use ark_ff::Zero;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use ark_std::test_rng;
use crate::sumcheck::{SumcheckProver, random_gkr_instance};
use ark_bls12_381::Fr;

pub struct NaiveSumCheck {
}

impl NaiveSumCheck {
    pub fn _new() -> NaiveSumCheck {
        return NaiveSumCheck {  };
    }

    pub fn convert_index_to_field<F: ark_ff::Field>(index: usize, field_size: usize) -> Vec<F> {
        todo!()
    }
}

impl<F: ark_ff::Field> SumcheckProver<F> for NaiveSumCheck {

    fn compute_sum(
            mult: &SparseMultilinearExtension<F>,
            vi: &DenseMultilinearExtension<F>,
            vj: &DenseMultilinearExtension<F>,
            r: &[F],
        ) {
        assert_eq!(mult.num_vars, vi.num_vars * 3);
        assert_eq!(mult.num_vars, vj.num_vars * 3);
        // num_vars is the same as the amount of bits needed to describe a gate.
        let mut sum= F::zero();
        
        // This fixes the gate label to the mult.
        let mult_gate = mult.fix_variables(r);
        
        // Can not use iterator as it iterates over the evaluations for inputs {0,1}^s.
        let total_bits= 1 << vj.num_vars; // Denotes the maxmimum bit value
        // Our for loop thus becomes for each i in {0,1}^s
        for i in 0..total_bits {
            let prefix_eval = vi[i];

            // We have to "bind" the next set of bits to the prefix in the mult.
            // First however we must convert "i" to a field so we can properly fix the "variables".
            let field_index = &NaiveSumCheck::convert_index_to_field(i, vj.num_vars);
            let mult_pref = mult_gate.fix_variables(field_index);

            for j in 0..total_bits {
                let suffix_eval = vj[j];
                sum += mult_pref[j] * prefix_eval * suffix_eval;
            }
        }
        todo!()
    }
}

#[test]
fn naive_sum_check_compute_sum() {
    let gate_labels_size = 2^4;
    let mut rng = test_rng();
    //let (mult, vi, vj) = random_gkr_instance::<ark_bls12_381::Fr>(gate_labels_size, &mut rng);
    assert_ne!(1, 0);
}