use std::time::Duration;
use ark_bls12_381::Fr;
use ark_std::test_rng;
use structures::data_structures::SumCheckProver;
use crate::gkr::gkr_round::GKRRound;
use crate::provers::fast::FastProver;
use crate::provers::naive::NaiveProver;
use crate::structures::circuit_structures::GKRCircuit;

mod util;
pub mod provers;
pub mod structures;
pub mod verifiers;
pub mod gkr;

fn main() {
    println!("Hey you just ran the main function!");
}

