# SumCheckSimulation

Simulation and benchmarking of the sum-check protocol — a naive prover
compared against a linear-time (prefix/suffix) prover — and of the GKR
protocol over layered arithmetic circuits. Built on arkworks (BLS12-381
scalar field).

## Running

The benchmark is selected with a command-line argument (no need to edit
`main.rs`):

```sh
# Sum-check benchmark with the config used in the report.
# Writes sum-check.xlsx. Note: took roughly 5 hours in total on my machine.
cargo run --release -- sumcheck

# GKR benchmark on the large random circuit used in the report.
# Writes gkr_circuit.xlsx.
cargo run --release -- gkr

# GKR benchmark on the small hand-built circuit from the report's figure.
cargo run --release -- gkr --small
```

## Project layout

```
src/
├─ main.rs        thin CLI entry point
├─ sumcheck/      the sum-check protocol
│  ├─ protocol.rs    SumCheckProver / SumCheckVerifier traits,
│  │                 layer-reduction message, line restriction
│  ├─ instance.rs    GKRRound: the per-layer instance a sum-check runs over
│  ├─ naive.rs       NaiveProver (recomputes the sum from scratch each round)
│  ├─ fast.rs        FastProver (linear time, folding bookkeeping tables)
│  ├─ verifier.rs    StandardVerifier
│  └─ tests.rs       shared prover tests, written once and instantiated
│                    for both provers
├─ gkr/           the GKR protocol
│  ├─ circuit.rs     Gate / Layer / GKRCircuit + evaluation, InputLayer
│  ├─ predicates.rs  per-layer ~add / ~mult wiring predicates
│  ├─ prover.rs      GKRProver
│  ├─ verifier.rs    GKRVerifier
│  └─ driver.rs      GKRDriver: runs one sum-check per layer + final check
├─ timing.rs      Track / AnalysisResult / timed() instrumentation
├─ bench/         benchmark harnesses + xlsx output
└─ util.rs        small shared helpers
```

## Output files

- `sumcheck` writes `sum-check.xlsx` (per-dimension prover/verifier times).
- `gkr` writes `gkr_circuit.xlsx` (per-layer prover/verifier times).

## Tests

```sh
cargo test --release
```

All randomness in tests and simulations uses the deterministically seeded
`ark_std::test_rng`, so nothing here is cryptographically secure — it is a
simulation for benchmarking, not a production proof system.
