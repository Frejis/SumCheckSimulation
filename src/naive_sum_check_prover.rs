use crate::{ark_sumcheck, rng::{Blake2b512Rng, FeedableRNG}, sumcheck::{random_gkr_instance, GKRRoundSumcheck, SumcheckProver}};
use ark_ff::{Field};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};
use ark_std::test_rng;

pub struct NaiveSumCheck {
}

impl NaiveSumCheck {
    pub fn new() -> NaiveSumCheck {
        return NaiveSumCheck {  };
    }

    // Converts an {0,1}^s for field size $s$ into a field element where each bit becomes a field.
    pub fn convert_index_to_field<F: Field>(index: usize, mut field_size: usize) -> Vec<F> {
        let mut result: Vec<F> = Vec::with_capacity(field_size);
        let mut tmp = index.clone();
        while field_size != 0 {
            let bit_value = (tmp & 1) as u64;
            result.push(bit_value.into());
            tmp = tmp >> 1;
            field_size -= 1;
        }
        result
    }
}

impl<F: Field> SumcheckProver<F> for NaiveSumCheck {

    fn compute_sum(
            mult: &SparseMultilinearExtension<F>,
            vi: &DenseMultilinearExtension<F>,
            vj: &DenseMultilinearExtension<F>,
            r: &[F],
        ) -> F {
        assert_eq!(mult.num_vars, vi.num_vars * 3);
        assert_eq!(mult.num_vars, vj.num_vars * 3);
        // num_vars is the same as the amount of bits needed to describe a gate.
        let mut sum= F::zero();
        let total_bits= 1 << vj.num_vars; // Denotes the maxmimum bit value

        // This fixes the gate label to the mult.
        let mult_gate = mult.fix_variables(r);
        
        // Can not use iterator as it iterates over the evaluations for inputs {0,1}^s.
        // Our for loop thus becomes for each i in {0,1}^s
        for i in 0..total_bits {
            let prefix_eval = vi[i];

            // We have to "bind" the next set of bits to the prefix in the mult.
            // First however we must convert "i" to a field so we can properly fix the "variables".
            //let field_index = &NaiveSumCheck::convert_index_to_field(i, vj.num_vars);
            let field_index = &NaiveSumCheck::convert_index_to_field(i, vj.num_vars);
            let mult_pref = mult_gate.fix_variables(field_index);

            for j in 0..total_bits {
                let suffix_eval = vj[j];
                sum += mult_pref[j] * prefix_eval * suffix_eval;
            }
        }
        sum
    }
}

mod tests {

#[test]
    fn test_my_naive_equals_arkworks() {
        super::test_my_naive_equal_ark_naive::<ark_bls12_381::Fr>(9);
    }
}

fn test_my_naive_equal_ark_naive<F: Field>(nv: usize) {
    let mut rng = test_rng();
    let (f1, f2, f3) = random_gkr_instance(nv, &mut rng);
    let g: Vec<_> = (0..nv).map(|_| F::rand(&mut rng)).collect();
    
    let my_claimed_sum: F = NaiveSumCheck::compute_sum(&f1, &f2, &f3, &g);
    let ark_claimed_sum: F = ark_sumcheck::calculate_sum_naive(&f1, &f2, &f3, &g);
    assert_eq!(my_claimed_sum, ark_claimed_sum);
}

fn test_naive<F: Field>(nv: usize) {
    let mut rng = test_rng();
    let (f1, f2, f3) = random_gkr_instance(nv, &mut rng);
    let g: Vec<_> = (0..nv).map(|_| F::rand(&mut rng)).collect();
    
    let claimed_sum = ark_sumcheck::calculate_sum_naive(&f1, &f2, &f3, &g);
    
    let mut rng = Blake2b512Rng::setup();
    let proof = GKRRoundSumcheck::prove(&mut rng, &f1, &f2, &f3, &g);
    rng = Blake2b512Rng::setup();
    let subclaim = GKRRoundSumcheck::verify(&mut rng, f2.num_vars, &proof, claimed_sum)
        .expect("verification failed");
    let result = subclaim.verify_subclaim(&f1, &f2, &f3, &g);
    assert!(result)
}