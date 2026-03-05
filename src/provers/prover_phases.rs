#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProverPhase {
    Uninitialized,
    PhaseOne,
    PhaseTwo,
}