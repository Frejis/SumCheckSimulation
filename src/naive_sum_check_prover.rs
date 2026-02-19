use std::ops::AddAssign;

use ark_ff::{Field, UniformRand, Zero};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial, SparseMultilinearExtension};
use ark_std::{rand::{Rng, RngCore}, test_rng};
use crate::{ark_sumcheck, rng::{Blake2b512Rng, FeedableRNG}, sumcheck::{GKRRoundSumcheck, SumcheckProver, random_gkr_instance}};
use ark_bls12_381::Fr;

pub struct NaiveSumCheck {
}

impl NaiveSumCheck {
    pub fn new() -> NaiveSumCheck {
        return NaiveSumCheck {  };
    }

    fn index_to_field_element<F: Field>(mut index: usize, mut nv: usize) -> Vec<F> {
        let mut ans = Vec::with_capacity(nv);
        while nv != 0 {
            ans.push(((index & 1) as u64).into());
            index >>= 1;
            nv -= 1;
        }
        ans
    }

/*
    pub fn convert_index_to_field<F: ark_ff::Field>(index: usize, field_size: usize) -> Vec<F> {
        let mut nv: Vec<F> = Vec::with_capacity(field_size);
        let mut tmp = index.clone();
        while index != 0 {
            let bit_value = (index & 1) as u64;
            nv.push(bit_value.into());
            tmp = tmp >> 1;
        }
        todo!()
    }
*/
}

impl<F: ark_ff::Field> SumcheckProver<F> for NaiveSumCheck {

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
            let field_index = &NaiveSumCheck::index_to_field_element(i, vj.num_vars);
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
    fn test_ark_naive() {
        super::test_my_naive::<ark_bls12_381::Fr>(9);

    }

#[test]
    fn test_my_naive() {
        super::test_my_naive::<ark_bls12_381::Fr>(9);
    }
}

fn test_my_naive<F: ark_ff::Field>(nv: usize) {
    let mut rng = test_rng();
    let (f1, f2, f3) = random_gkr_instance(nv, &mut rng);
    let g: Vec<_> = (0..nv).map(|_| F::rand(&mut rng)).collect();
    
    let claimed_sum = NaiveSumCheck::compute_sum(&f1, &f2, &f3, &g);
    
    let mut rng = Blake2b512Rng::setup();
    let proof = GKRRoundSumcheck::prove(&mut rng, &f1, &f2, &f3, &g);
    rng = Blake2b512Rng::setup();
    let subclaim = GKRRoundSumcheck::verify(&mut rng, f2.num_vars, &proof, claimed_sum)
        .expect("verification failed");
    let result = subclaim.verify_subclaim(&f1, &f2, &f3, &g);
    assert!(result)
}

fn test_naive<F: ark_ff::Field>(nv: usize) {
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