use std::fmt::format;
use std::time::{Duration, Instant};
use ark_bls12_381::{Fr, FrConfig};
use ark_ff::{Field, Fp, MontBackend, Zero};
use ark_poly::univariate::SparsePolynomial;
use ark_std::test_rng;
use rust_xlsxwriter::{ColNum, RowNum, Workbook};
use crate::gkr::gkr_circuit::compute_predicates;
use crate::gkr::gkr_driver::GKRDriver;
use crate::gkr::gkr_prover::GKRProver;
use crate::gkr::gkr_round::GKRRound;
use crate::gkr::gkr_verifier::GKRVerifier;
use crate::gkr::layer::InputLayer;
use crate::provers::fast::FastProver;
use crate::provers::naive::NaiveProver;
use crate::structures::circuit_structures::GKRCircuit;
use crate::structures::data_structures::{AnalysisResult, SumCheckProver, SumCheckVerifier, Track};
use crate::util::{create_prover, random_gate, sparse_polynomial};
use crate::verifiers::standard_verifier::StandardVerifier;

mod util;
pub mod provers;
pub mod structures;
pub mod verifiers;
pub mod gkr;

fn main() {
    benchmark_gkr();
}

fn run_sum_check_config() {
    let mut trials: Vec<usize> = vec![0; 6];
    trials.append(&mut vec![300; 5]);
    trials.append(&mut vec![30; 10]);
    simulate_circuit_and_save_results(14, &*trials)
}

fn simulate_circuit_and_save_results(max_variables: usize, trials: &[usize]) {
    let mut naive_res: Vec<Vec<Track>> = Vec::new();
    let mut fast_res: Vec<Vec<Track>> = Vec::new();
    for i in 0..max_variables {
        let trial = trials[i];
        let mut fast_trials = Vec::new();
        let mut naive_trials = Vec::new();
        for _ in 0..trial {
            let naive_time = simulate_sum_check_instance::<NaiveProver<Fr>, Fr, StandardVerifier<Fr>>(i);
            println!("Finished running naive for dimension {i}");
            println!("Time taken for naive prover: {:?}", naive_time.prover());
            naive_trials.push(naive_time);

            let fast_time = simulate_sum_check_instance::<FastProver<Fr>, Fr, StandardVerifier<Fr>>(i);
            println!("Finished running fast for dimension {i}");
            println!("Time taken for fast prover: {:?}", fast_time.prover());
            fast_trials.push (fast_time);
        }
        naive_res.push(naive_trials);
        fast_res.push(fast_trials);
    }
    let mut workbook = Workbook::new();
    let worksheet = workbook.add_worksheet();
    // Create the headers.
    for i in 0..max_variables {
        worksheet.write(0, (i * 4) as ColNum, format!("Fast Prover, {i} variables")).unwrap();
        worksheet.write(0, (i * 4 + 1) as ColNum, format!("Fast Verifier, {i} variables")).unwrap();
        worksheet.write(0, (i * 4 + 2) as ColNum, format!("Naive Prover, {i} variables")).unwrap();
        worksheet.write(0, (i * 4 + 3) as ColNum, format!("Naive Verifier, {i} variables")).unwrap();
    }

    for i in 0..naive_res.len() {
        for j in 0..naive_res[i].len() {
            let naive_prover_time = naive_res[i][j].prover().as_secs_f64();
            let naive_verifier_time = naive_res[i][j].verifier().as_secs_f64();
            let fast_prover_time = fast_res[i][j].prover().as_secs_f64();
            let fast_verifier_time = fast_res[i][j].verifier().as_secs_f64();

            let row = (1 + j) as RowNum;

            worksheet.write(row, (i * 4) as ColNum, fast_prover_time).unwrap();
            worksheet.write(row, (i * 4 + 1) as ColNum, fast_verifier_time).unwrap();
            worksheet.write(row, (i * 4 + 2) as ColNum, naive_prover_time).unwrap();
            worksheet.write(row, (i * 4 + 3) as ColNum, naive_verifier_time).unwrap();
        }
    }

    workbook.save("sum-check.xlsx").unwrap();
    println!("Saved data to excel sheet");
}

fn simulate_sum_check_instance<T: SumCheckProver<F>, F:Field, V: SumCheckVerifier<F>>(dim: usize) -> Track {
    let r_gkr: GKRRound<F> = GKRRound::new_rand(dim);
    let r_gate = random_gate(r_gkr.gate_labes());
    let mut track = Track::new();

    let (mut prv, time) = create_prover::<T, F>(&mut r_gate.clone(), r_gkr.clone());
    track.add_prover_time(time);
    let (claim, time) = prover_compute_claim::<T, F>(&mut prv);
    track.add_prover_time(time);

    let (mut vrf, time) = create_verifier::<F, V>(claim);
    track.add_verifier_time(time);

    for i in 0..r_gkr.vj.num_vars*2 {
        let (vrf_fnc, time) = track_compute_vrf_func(&mut prv);
        track.add_prover_time(time);

        if i == r_gkr.vj.num_vars*2-1 {
            let tmp = vrf_fnc.iter().rfold(F::zero(), |acc, (_, elem)| {*elem + acc});
            assert_eq!(tmp, prv.compute_sum());
            let ins = Instant::now();
            //vrf.final_check(&*r_gate, r_gkr.clone().add_predicate(), r_gkr.clone().mult_predicate(), r_gkr.vi, vrf_func);
            track.add_verifier_time(ins.elapsed());
            break;
        }

        let (r_i, time) = handle_vrf_checks(&mut vrf, vrf_fnc.clone());
        track.add_verifier_time(time);

        let time = fix_variable_prv(r_i, &mut prv);
        track.add_prover_time(time);

        let (claim, time) = prover_compute_claim(&mut prv);
        track.add_prover_time(time);

        let time = set_new_claim(&mut vrf, claim);
        track.add_verifier_time(time);

    }
    track
}

fn set_new_claim<F: Field, V: SumCheckVerifier<F>>(vrf: &mut V, claim: F) -> Duration {
    let ins = Instant::now();
    vrf.set_claim(claim);
    ins.elapsed()
}

fn fix_variable_prv<T: SumCheckProver<F>, F: Field>(r_i: F, prv: &mut T) -> Duration {
    let ins = Instant::now();
    prv.fix_variable(r_i);
    ins.elapsed()
}

fn handle_vrf_checks<F: Field, V: SumCheckVerifier<F>>(vrf: &mut V, evals: Vec<(usize, F)>) -> (F, Duration) {
    let ins = Instant::now();
    let polynomial = SparsePolynomial::from_coefficients_vec(evals);
    let r_i = vrf.handle_round(&polynomial);
    (r_i, ins.elapsed())
}

fn track_compute_vrf_func<T: SumCheckProver<F>, F: Field>(prv: &mut T) -> (Vec<(usize, F)>, Duration) {
    let ins = Instant::now();
    let vrf_func = prv.get_verifier_function();
    (vrf_func, ins.elapsed())
}

fn create_verifier<F: Field, V: SumCheckVerifier<F>>(claim: F) -> (V, Duration) {
    let ins = Instant::now();
    // The degree is abitrary atm since i am unsure how to get the maximum degree for the random
    // polynomials sent during Sum-check.
    let vrf = V::new(2, claim);
    (vrf, ins.elapsed())
}

fn prover_compute_claim<T: SumCheckProver<F>, F: Field>(prv: &mut T) -> (F, Duration) {
    let inst = Instant::now();
    let claim = prv.compute_sum();
    (claim, inst.elapsed())
}

fn benchmark_gkr() {
    let (layers, random_circuit, input_layer) = random_circuit();
    let mut naive_res: Vec<AnalysisResult> = Vec::new();
    let mut fast_res: Vec<AnalysisResult> = Vec::new();
    let trials = 30;
    for _ in 0..trials {
        let naive = simulate_gkr_naive::<Fr>(random_circuit.clone(), input_layer.clone());
        let fast = simulate_gkr_fast::<Fr>(random_circuit.clone(), input_layer.clone());
        naive_res.push(naive);
        fast_res.push(fast)
    }

    let mut fast_time: Vec<Vec<Track>> = Vec::new();
    let mut naive_time: Vec<Vec<Track>> = Vec::new();
    for i in 0..layers.len() {
        let mut fast_for_layer: Vec<Track> = Vec::new();
        let mut naive_for_layer: Vec<Track> = Vec::new();
        for j in 0..naive_res.len() {
            fast_for_layer.push(fast_res[j].get_time_for_layer(i).clone());
            naive_for_layer.push(naive_res[j].get_time_for_layer(i).clone());
        }
        fast_time.push(fast_for_layer);
        naive_time.push(naive_for_layer);
    }
    save_results_for_each_layer_avg(fast_time, naive_time, layers);
}

pub fn save_results_for_each_layer_avg(fast: Vec<Vec<Track>>, naive: Vec<Vec<Track>>, layers: Vec<usize>) {
    let mut workbook = Workbook::new();
    let worksheet = workbook.add_worksheet();
    for i in 0..layers.len() {
        worksheet.write(0,(i * 5) as ColNum, "Variables/Layer").unwrap();
        worksheet.write(0, (i * 5 + 1) as ColNum, format!("Fast Prover layer {i}")).unwrap();
        worksheet.write(0, (i * 5 + 2) as ColNum, format!("Fast Verifier layer {i}")).unwrap();
        worksheet.write(0, (i * 5 + 3) as ColNum, format!("Naive Prover layer {i}")).unwrap();
        worksheet.write(0, (i * 5 + 4) as ColNum, format!("Naive Verifier layer {i}")).unwrap();
    }

    for i in 0..layers.len() {
        worksheet.write((i + 1) as RowNum, 0, layers[i].ilog2()).unwrap();
        for j in 0..fast[i].len() {
            let row = j + 1;
            worksheet.write(row as RowNum, (i * 5 + 1) as ColNum, fast[i][j].prover().as_secs_f64()).unwrap();
            worksheet.write(row as RowNum, (i * 5 + 2) as ColNum, fast[i][j].verifier().as_secs_f64()).unwrap();
            worksheet.write(row as RowNum, (i * 5 + 3) as ColNum, naive[i][j].prover().as_secs_f64()).unwrap();
            worksheet.write(row as RowNum, (i * 5 + 4) as ColNum, naive[i][j].verifier().as_secs_f64()).unwrap();
        }
    }
    workbook.save("gkr_circuit.xlsx").unwrap();
    println!("Saved data to excel sheet");
}

fn random_circuit() -> (Vec<usize>, GKRCircuit<Fr>, InputLayer<Fr>) {
    //let layers = &vec![2, 4, 8, 32, 64, 128, 256, 512, 2048, 2048*2, 2048*8, 1024];
    let layers = &vec![2048, 2048*2, 2048*8, 2048*8, 1024];
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

