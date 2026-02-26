// This is a prover. It takes a circuit and a sumcheck helper and proves the circuit for the verifier.

use ark_ff::Field;
use ark_poly::DenseMultilinearExtension;
use crate::circuit_structures::GkrCircuit;
use crate::data_structures::SumCheckProver;

struct Prover<F: Field, P: SumCheckProver<F>> {
    sc_prover: P,
    circuit: GkrCircuit<F>,
}

impl<F: Field, P: SumCheckProver<F>> Prover<F, P> {
    pub fn new(sc_prover: P, circuit: GkrCircuit<F>) -> Self {
        Self {
            sc_prover,
            circuit,
        }
    }

    pub fn sum(&mut self) -> F {
        self.sc_prover.compute_sum()
    }

    pub fn fix_variable_sum(&mut self, f: F) {
        self.sc_prover.fix_variable(f);
    }

    pub fn fix_variable_last_round(&mut self) -> DenseMultilinearExtension<F> {
        let z1 = self.sc_prover.compute_z_1();
        let z2 = self.sc_prover.compute_z_2();
        DenseMultilinearExtension::from_evaluations_vec(1, vec![z1, z2])
    }
}