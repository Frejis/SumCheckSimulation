use ark_bls12_381::Fr;
use ark_std::test_rng;
use crate::gkr::gkr_driver::GKRDriver;
use crate::gkr::gkr_prover::GKRProver;
use crate::gkr::gkr_verifier::GKRVerifier;
use crate::gkr::layer::InputLayer;
use crate::structures::circuit_structures::GKRCircuit;

mod util;
pub mod provers;
pub mod structures;
pub mod verifiers;
pub mod gkr;

fn main() {
    let layers = &[2, 4, 8, 32, 64];
    let random_circuit: GKRCircuit<Fr> = GKRCircuit::random(layers, &mut test_rng());
    let input_layer: InputLayer<Fr> = InputLayer::random(layers.last().unwrap());
    let mut gkr_prover = GKRProver::new(random_circuit.clone(), input_layer.clone());
    gkr_prover.compute_predicates();
    let gkr_verifier = GKRVerifier::new(random_circuit.clone(), input_layer.clone());
    let mut gkrdriver: GKRDriver<Fr> = GKRDriver::new(gkr_prover, gkr_verifier, random_circuit, input_layer);
    gkrdriver.run_circuit();
    println!("Process did not crash");
}

