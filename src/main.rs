use std::time::Duration;
use ark_bls12_381::{Fr, FrConfig};
use ark_ff::{Field, Fp, MontBackend};
use ark_std::test_rng;
use rust_xlsxwriter::{RowNum, Workbook};
use crate::gkr::gkr_circuit::compute_predicates;
use crate::gkr::gkr_driver::GKRDriver;
use crate::gkr::gkr_prover::GKRProver;
use crate::gkr::gkr_verifier::GKRVerifier;
use crate::gkr::layer::InputLayer;
use crate::provers::fast::FastProver;
use crate::provers::naive::NaiveProver;
use crate::structures::circuit_structures::GKRCircuit;
use crate::structures::data_structures::{AnalysisResult, SumCheckProver, Track};

mod util;
pub mod provers;
pub mod structures;
pub mod verifiers;
pub mod gkr;

fn main() {
    let (layers, random_circuit, input_layer) = random_circuit();
    let mut naive_res: Vec<AnalysisResult> = Vec::new();
    let mut fast_res: Vec<AnalysisResult> = Vec::new();
    let trials = 8;
    for _ in 0..trials {
        let naive = simulate_gkr_naive::<Fr>(random_circuit.clone(), input_layer.clone());
        let fast = simulate_gkr_fast::<Fr>(random_circuit.clone(), input_layer.clone());
        naive_res.push(naive);
        fast_res.push(fast)
    }

    let avg_pr_layer_naive = compute_avg_layers(layers.clone(), &mut naive_res, trials);
    let avg_pr_layer_fast = compute_avg_layers(layers.clone(), &mut fast_res, trials);
    println!("Naive prover:");
    for (layer, time) in avg_pr_layer_naive.iter().enumerate() {
        println!("Average Time for layer {layer}. Prover {:?}, Verifier {:?}", time.prover(), time.verifier());
    }
    println!("Fast prover:");
    for (layer, time) in avg_pr_layer_fast.iter().enumerate() {
        let variables = layers[layer].ilog2();
        println!("Average Time for layer {layer}. Prover {:?}, Verifier {:?}. Variables in layer: {:?}", time.prover(), time.verifier(), variables);
    }

    save_results(avg_pr_layer_fast, avg_pr_layer_naive, layers);
}

pub fn save_results(fast: Vec<Track>, naive: Vec<Track>, layers: Vec<usize>) {
    let mut workbook = Workbook::new();
    let worksheet = workbook.add_worksheet();
    worksheet.write(0,0, "Variables/Layer").unwrap();
    worksheet.write(0, 1, "Fast Prover").unwrap();
    worksheet.write(0, 2, "Fast Verifier").unwrap();
    worksheet.write(0, 3, "Naive Prover").unwrap();
    worksheet.write(0, 4, "Naive Verifier").unwrap();

    for i in 0..fast.len() {
        let row = i + 1;
        worksheet.write(row as RowNum, 0, layers[i].ilog2()).unwrap();
        worksheet.write(row as RowNum, 1, fast[i].prover().as_secs_f64()).unwrap();
        worksheet.write(row as RowNum, 2, fast[i].verifier().as_secs_f64()).unwrap();
    }
    for i in 0..naive.len() {
        let row = i + 1;
        worksheet.write(row as RowNum, 3, naive[i].prover().as_secs_f64()).unwrap();
        worksheet.write(row as RowNum, 4, naive[i].verifier().as_secs_f64()).unwrap();
    }
    workbook.save("test.xlsx").unwrap();
    println!("Saved data to excel sheet");
}

fn random_circuit() -> (Vec<usize>, GKRCircuit<Fr>, InputLayer<Fr>) {
    //let layers = &vec![2, 4, 8, 32, 64, 128, 256, 512, 2048, 2048*2, 2048*8, 1024];
    let layers = &vec![2048, 2048*2, 2048*8, 1024];
    let random_circuit: GKRCircuit<Fr> = GKRCircuit::random(layers, &mut test_rng());
    let input_layer: InputLayer<Fr> = InputLayer::random(layers.last().unwrap());
    (layers.clone(), random_circuit, input_layer)
}

fn figure_circuit() -> (Vec<usize>, GKRCircuit<Fr>, InputLayer<Fr>) {
    let layers = &vec![2, 4, 4];
    let circuit: GKRCircuit<Fr> = GKRCircuit::figure_circuit();
    let input_layer: InputLayer<Fr> = InputLayer::random(layers.last().unwrap());
    (layers.clone(), circuit, input_layer)
}

fn compute_avg_layers(layers: Vec<usize>, analysis: &mut Vec<AnalysisResult>, trials: usize) -> Vec<Track> {
    let mut track = Vec::new();
    for i in 0..layers.len() {
        let mut prv_time = Duration::new(0, 0);
        let mut vrf_time = Duration::new(0, 0);

        for j in 0..trials {
            prv_time += analysis[j].get_time_for_layer(i).prover();
            vrf_time += analysis[j].get_time_for_layer(i).verifier();
        }
        prv_time = prv_time / trials as u32;
        vrf_time = vrf_time / trials as u32;
        track.push(Track::new_times(prv_time, vrf_time));
    }
    track
}

fn simulate_gkr_naive<F: Field>(random_circuit: GKRCircuit<Fr>, input_layer: InputLayer<Fr>)
                                -> AnalysisResult
where NaiveProver<F>: SumCheckProver<Fp<MontBackend<FrConfig, 4>, 4>> {
    let mut gkr_prover = GKRProver::new(random_circuit.clone(), input_layer.clone());
    let predicate = compute_predicates(gkr_prover.circuit(), gkr_prover.input());
    gkr_prover.set_predicates(predicate.clone());
    let mut gkr_verifier = GKRVerifier::new(random_circuit.clone(), input_layer.clone());
    gkr_verifier.set_predicate(predicate);
    let mut gkrdriver: GKRDriver<Fr> = GKRDriver::new(&gkr_prover, &gkr_verifier, random_circuit, input_layer);
    gkrdriver.run_circuit::<NaiveProver<F>>()
}

fn simulate_gkr_fast<F: Field>(random_circuit: GKRCircuit<Fr>, input_layer: InputLayer<Fr>)
    -> AnalysisResult
where FastProver<F>: SumCheckProver<Fp<MontBackend<FrConfig, 4>, 4>> {
    let mut gkr_prover = GKRProver::new(random_circuit.clone(), input_layer.clone());
    let mut gkr_verifier = GKRVerifier::new(random_circuit.clone(), input_layer.clone());
    let predicate = compute_predicates(gkr_prover.circuit(), gkr_prover.input());
    gkr_prover.set_predicates(predicate.clone());
    gkr_verifier.set_predicate(predicate);
    
    
    let mut gkrdriver: GKRDriver<Fr> = GKRDriver::new(&gkr_prover, &gkr_verifier, random_circuit, input_layer);
    gkrdriver.run_circuit::<FastProver<F>>()
}

