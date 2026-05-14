use std::time::Duration;
use ark_bls12_381::{Fr, FrConfig};
use ark_ff::{Field, Fp, MontBackend};
use ark_std::test_rng;
use crate::gkr::gkr_driver::GKRDriver;
use crate::gkr::gkr_prover::GKRProver;
use crate::gkr::gkr_verifier::GKRVerifier;
use crate::gkr::layer::InputLayer;
use crate::provers::fast::FastProver;
use crate::provers::naive::NaiveProver;
use crate::structures::circuit_structures::GKRCircuit;
use crate::structures::data_structures::{AnalysisResult, SumCheckProver};

mod util;
pub mod provers;
pub mod structures;
pub mod verifiers;
pub mod gkr;

fn main() {
    let layers = &[2, 4, 8, 32, 64, 128, 256, 512, 2048, 1024];
    let random_circuit: GKRCircuit<Fr> = GKRCircuit::random(layers, &mut test_rng());
    let input_layer: InputLayer<Fr> = InputLayer::random(layers.last().unwrap());
    let naive_res = simulate_gkr_naive::<Fr>(random_circuit.clone(), input_layer.clone());
    let fast_res = simulate_gkr_fast::<Fr>(random_circuit, input_layer);
    println!("Process did not crash");
    println!("Naive results:");
    naive_res.print();
    println!("Fast results:");
    fast_res.print();
}

fn simulate_gkr_naive<F: Field>(random_circuit: GKRCircuit<Fr>, input_layer: InputLayer<Fr>)
    -> AnalysisResult
where NaiveProver<F>: SumCheckProver<Fp<MontBackend<FrConfig, 4>, 4>> {
    let mut gkr_prover = GKRProver::new(random_circuit.clone(), input_layer.clone());
    gkr_prover.compute_predicates();
    let gkr_verifier = GKRVerifier::new(random_circuit.clone(), input_layer.clone());
    let mut gkrdriver: GKRDriver<Fr> = GKRDriver::new(gkr_prover, gkr_verifier, random_circuit, input_layer);
    gkrdriver.run_circuit::<NaiveProver<F>>()
}

fn simulate_gkr_fast<F: Field>(random_circuit: GKRCircuit<Fr>, input_layer: InputLayer<Fr>)
    -> AnalysisResult
where FastProver<F>: SumCheckProver<Fp<MontBackend<FrConfig, 4>, 4>> {
    let mut gkr_prover = GKRProver::new(random_circuit.clone(), input_layer.clone());
    gkr_prover.compute_predicates();
    let gkr_verifier = GKRVerifier::new(random_circuit.clone(), input_layer.clone());
    let mut gkrdriver: GKRDriver<Fr> = GKRDriver::new(gkr_prover, gkr_verifier, random_circuit, input_layer);
    gkrdriver.run_circuit::<FastProver<F>>()
}

