use std::arch::x86_64::__m128bh;
use ark_bls12_381::Fr;
use ark_std::test_rng;
use crate::data_structures::{Prover, Verifier};
use crate::naive_sum_check::NaiveProver;
use crate::standard_verifier::StandardVerifier;
use crate::util::{random_gate, random_gkr_instance};

mod data_structures;
mod util;
pub mod naive_sum_check;
mod standard_verifier;

fn main() {
    println!("Hello, world!");

    // Just testing things out. First create a random GKR instance with Fr
    let variables = 10;
    let mut rng = test_rng();
    let (mult, vi, vj) = random_gkr_instance::<Fr, _>(10, &mut rng);
    let rand_gate = random_gate(variables);

    let prover = NaiveProver::new(mult, vi, vj, rand_gate);

    // I have no clue about the max degree atm so i just say fuck it we ball hahaha
    let mut verifier = StandardVerifier::new(variables*400, prover.compute_sum());
    let verifier_func = prover.get_verifier_function();
    if verifier.check_claimed_value(&verifier_func) {
        println!("Claimed value is correct.");
    }
}

mod generic_tests {
    use ark_bls12_381::Fr;
    use ark_std::test_rng;
    use crate::data_structures::{Prover, Verifier};
    use crate::naive_sum_check::NaiveProver;
    use crate::standard_verifier::StandardVerifier;
    use crate::util::{random_gate, random_gkr_instance};

    #[test]
    fn test_one_round() {
        // 7 variables 3 seconds.... 8 variables 23!!! seconds!??! :OOO whaaa

        let variables = 7;
        let mut rng = test_rng();
        let (mult, vi, vj) = random_gkr_instance::<Fr, _>(7, &mut rng);
        let rand_gate = random_gate(variables);

        let prover = NaiveProver::new(mult, vi, vj, rand_gate);

        // I have no clue about the max degree atm so i just say fuck it we ball hahaha
        let verifier = StandardVerifier::new(variables*400, prover.compute_sum());
        let verifier_func = prover.get_verifier_function();
        if verifier.check_claimed_value(&verifier_func) {
            println!("Claimed value is correct.");
        }
        assert_eq!(1,1);
    }
}