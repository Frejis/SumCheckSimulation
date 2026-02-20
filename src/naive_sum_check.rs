use ark_ff::{Field, Zero};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};
use rand::random;
use crate::data_structures::Prover;
use crate::util::index_to_field_element;

struct NaiveProver<F: Field> {
    mult: SparseMultilinearExtension<F>, // mle of mult for round k with (r, i , j)
    vi: DenseMultilinearExtension<F>, // mle of v_(k-1)(i)
    vj: DenseMultilinearExtension<F>, // mle of v_(k-1)(i)
    r: Vec<F>, // The gate "r" that is fixed.
}

impl<F: Field> NaiveProver<F> {
    pub fn new(
        mult: SparseMultilinearExtension<F>,
        vi: DenseMultilinearExtension<F>,
        vj: DenseMultilinearExtension<F>,
        r: Vec<F>,
    ) -> NaiveProver<F> {
        assert_eq!(mult.num_vars, vi.num_vars * 3);
        assert_eq!(mult.num_vars, vj.num_vars * 3);
        NaiveProver {
            mult,
            vi,
            vj,
            r,
        }
    }
}

impl<F: Field> Prover<F> for NaiveProver<F> {
    // Needs to be refactored just my last sumcheck which i know works.
    fn compute_sum(&self) -> F {
        let mult = &self.mult;
        let vi = &self.vi;
        let vj = &self.vj;
        let r = &self.r;
        // num_vars is the same as the amount of bits needed to describe a gate.
        let mut sum= F::zero();
        let total_bits= 1 << vj.num_vars; // Denotes the maximum bit value

        // This fixes the gate label to the mult.
        let mult_gate = mult.fix_variables(&*r);

        // Can not use iterator as it iterates over the evaluations for inputs {0,1}^s.
        // Our for loop thus becomes for each i in {0,1}^s
        for i in 0..total_bits {
            let prefix_eval = vi[i];

            // We have to "bind" the next set of bits to the prefix in the mult.
            // First however we must convert "i" to a field so we can properly fix the "variables".
            //let field_index = &NaiveSumCheck::convert_index_to_field(i, vj.num_vars);
            let field_index = index_to_field_element(i, vj.num_vars);
            let mult_pref = mult_gate.fix_variables(&*field_index);

            for j in 0..total_bits {
                let suffix_eval = vj[j];
                sum += mult_pref[j] * prefix_eval * suffix_eval;
            }
        }
        sum
    }

    fn get_verifier_function(&self) -> DenseMultilinearExtension<F> {
        // clone existing functions.
        /*
        Iterate over all possible assignments of the bits.
        After that take the first variable and set it to 0 and the other one 1.
        */
        // Assume that the gate has been fixed.
        let n = self.vi.num_vars + self.vj.num_vars;
        assert_eq!(self.mult.num_vars, self.vi.num_vars + self.vj.num_vars);
        let remaining_variables = n - 1;
        let total = 1usize.checked_shl(remaining_variables as u32).expect("too many vars");
        let mut s0 = F::zero();
        let mut s1 = F::zero();

        for mask in 0..total {
            let mut field_index: Vec<F> = index_to_field_element(mask, n);

            // Set the first variable to 0.
            field_index[0] = F::zero();
            // Now evaluate the function and add it to s0.
            let v1 = eval_g(&field_index);
            s0 += &v1;

            field_index[0] = F::one();
            let v2 = eval_g(&field_index);
            s1 += &v2;
        }

        DenseMultilinearExtension::from_evaluations_vec(1, vec![s0, s1])
    }

    fn fix_variable(&mut self, random_field_element: F) {
        /*
        1. Fix the first variable in mult. Then fix in vi. Once vi has no more variables fix vj.
        */
        let field_packed = &[random_field_element];
        self.mult.fix_variables(field_packed);

        let vi = &self.vi;
        if (vi.num_vars > 0) {
            vi.fix_variables(field_packed);
        } else {
            self.vj.fix_variables(field_packed);
        }
    }
}

fn eval_g<F: Field>(p0: &Vec<F>) -> F {
    todo!()
}

mod tests {
    use ark_ff::{One, Zero};
    use crate::util::index_to_field_element;

    #[test]
    fn sanity_check() {
        for mask in 0..30 {
            let mut points = vec![ark_bls12_381::Fr::zero(); 30];
            let field_index: Vec<ark_bls12_381::Fr> = index_to_field_element(mask, 30);
            for j in 0..30 {
                let bit = (mask >> j) & 1 != 0;
                points[j] = if bit { ark_bls12_381::Fr::one() } else { ark_bls12_381::Fr::zero() }
            }
            assert_eq!(field_index, points)
        }
    }
}