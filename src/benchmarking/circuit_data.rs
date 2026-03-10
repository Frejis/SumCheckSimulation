use std::fs;
use std::path::{Path, PathBuf};
use ark_ff::Field;
use serde::{Deserialize, Serialize};
use crate::structures::circuit_structures::Gate;

/// Serializable representation of a layer
#[derive(Serialize, Deserialize)]
pub struct SerializableLayer {
    pub gates: Vec<Gate>,
    pub values: Vec<Vec<u8>>, // Field elements serialized to bytes
}

/// Serializable representation of a circuit
#[derive(Serialize, Deserialize)]
pub struct SerializableCircuit {
    pub layers: Vec<SerializableLayer>,
}

impl SerializableCircuit {
    /// Convert a GKRCircuit to a serializable format
    pub fn from_circuit<F: Field>(circuit: &crate::structures::circuit_structures::GKRCircuit<F>) -> Self {
        let layers = circuit
            .layers
            .iter()
            .map(|layer| {
                let gates = layer.gates.clone();
                let values = layer
                    .values
                    .iter()
                    .map(|val| {
                        let mut bytes = Vec::new();
                        val.serialize_compressed(&mut bytes).expect("Serialization failed");
                        bytes
                    })
                    .collect();
                SerializableLayer { gates, values }
            })
            .collect();
        SerializableCircuit { layers }
    }

    /// Convert back to a GKRCircuit
    pub fn to_circuit<F: Field>(&self) -> crate::structures::circuit_structures::GKRCircuit<F> {
        use crate::gkr::layer::Layer;
        use crate::structures::circuit_structures::GKRCircuit;

        let layers = self
            .layers
            .iter()
            .map(|ser_layer| {
                let gates = ser_layer.gates.clone();
                let values = ser_layer
                    .values
                    .iter()
                    .map(|bytes| {
                        F::deserialize_compressed(&bytes[..]).expect("Deserialization failed")
                    })
                    .collect();
                Layer { gates, values }
            })
            .collect();
        GKRCircuit { layers }
    }
}

/// Manages a collection of pre-generated circuit instances for benchmarking
pub struct BenchmarkCircuitSet {
    circuits: Vec<SerializableCircuit>,
}

impl BenchmarkCircuitSet {
    /// Create a new empty set
    pub fn new() -> Self {
        Self {
            circuits: Vec::new(),
        }
    }

    /// Add a circuit to the set
    pub fn add_circuit<F: Field>(&mut self, circuit: &crate::structures::circuit_structures::GKRCircuit<F>) {
        self.circuits.push(SerializableCircuit::from_circuit(circuit));
    }

    /// Get the number of circuits in the set
    pub fn len(&self) -> usize {
        self.circuits.len()
    }

    /// Check if the set is empty
    pub fn is_empty(&self) -> bool {
        self.circuits.is_empty()
    }

    /// Get a circuit by index
    pub fn get_circuit<F: Field>(&self, index: usize) -> Option<crate::structures::circuit_structures::GKRCircuit<F>> {
        self.circuits.get(index).map(|ser| ser.to_circuit())
    }

    /// Save the circuit set to a file
    pub fn save_to_file<P: AsRef<Path>>(&self, path: P) -> Result<(), String> {
        let encoded = bincode::serialize(&self.circuits)
            .map_err(|e| format!("Failed to serialize circuits: {}", e))?;

        // Create parent directory if it doesn't exist
        if let Some(parent) = path.as_ref().parent() {
            fs::create_dir_all(parent)
                .map_err(|e| format!("Failed to create directory: {}", e))?;
        }

        fs::write(path, encoded)
            .map_err(|e| format!("Failed to write file: {}", e))?;
        Ok(())
    }

    /// Load a circuit set from a file
    pub fn load_from_file<P: AsRef<Path>>(path: P) -> Result<Self, String> {
        let data = fs::read(path)
            .map_err(|e| format!("Failed to read file: {}", e))?;

        let circuits: Vec<SerializableCircuit> = bincode::deserialize(&data)
            .map_err(|e| format!("Failed to deserialize circuits: {}", e))?;

        Ok(Self { circuits })
    }

    /// Get an iterator over all circuits
    pub fn iter<F: Field>(&self) -> impl Iterator<Item = crate::structures::circuit_structures::GKRCircuit<F>> + '_ {
        self.circuits.iter().map(|ser| ser.to_circuit())
    }
}

impl Default for BenchmarkCircuitSet {
    fn default() -> Self {
        Self::new()
    }
}

/// Generate a default path for benchmark data
pub fn default_benchmark_data_path(layer_config: &[usize], num_trials: usize) -> PathBuf {
    let layer_str = layer_config
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>()
        .join("_");

    PathBuf::from("benchmark_data")
        .join(format!("circuits_{}_{}.bin", layer_str, num_trials))
}

