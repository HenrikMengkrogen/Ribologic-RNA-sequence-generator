use std::path::{Path, PathBuf};
use std::process::Command;

use crate::config::USE_ML_MODEL;

pub fn call_model(input_file: &str) -> Result<(), Box<dyn std::error::Error>> {
    if USE_ML_MODEL {
        println!("Calling ML model to generate seed sequences...");
    } else {
        println!("ML model usage is disabled. Skipping seed generation.");
        return Ok(());
    }

    let project_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    let python = project_dir.join(".venv").join("bin").join("python");
    let script = project_dir.join("Python").join("RUN_MODEL.py");
    let input = project_dir.join(input_file);
    let output = project_dir.join("Python").join("seeds.txt");

    eprintln!("Project directory: {}", project_dir.display());
    eprintln!("Python executable: {}", python.display());
    eprintln!("Python exists:     {}", python.exists());
    eprintln!("Model script:      {}", script.display());
    println!("Script exists:     {}", script.exists());
    println!("Input file:        {}", input.display());
    println!("Input exists:      {}", input.exists());
    println!("Output file:       {}", output.display());

    if !python.exists() {
        return Err(format!("Python virtual environment not found: {}", python.display()).into());
    }

    if !script.exists() {
        return Err(format!("RUN_MODEL.py not found: {}", script.display()).into());
    }

    if !input.exists() {
        return Err(format!("input.txt not found: {}", input.display()).into());
    }

    let status = Command::new(&python)
        .arg(&script)
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .arg("--n")
        .arg("16")
        .current_dir(&project_dir)
        .status()
        .map_err(|err| format!("Could not start Python at {}: {err}", python.display()))?;

    match status.code() {
        Some(0) => Ok(()),

        Some(1) => {
            eprintln!("Warning: model skipped one or more targets.");
            Ok(())
        }

        Some(2) => Err("ML model failed fatally; no usable seeds were generated.".into()),

        code => Err(format!("RUN_MODEL.py exited unexpectedly: {code:?}").into()),
    }
}

pub fn import_seed_sequences(seed_file: &Path) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    let contents = std::fs::read_to_string(seed_file)?;

    let seeds: Vec<String> = contents
        .lines()
        .filter_map(|line| {
            let (key, value) = line.split_once(':')?;

            if key.trim().eq_ignore_ascii_case("sequence") {
                let sequence = value.trim();

                if !sequence.is_empty() {
                    return Some(sequence.to_string());
                }
            }

            None
        })
        .collect();

    if seeds.is_empty() {
        return Err(format!("No seed sequences found in {}", seed_file.display()).into());
    }

    Ok(seeds)
}

pub fn train_model() -> Result<(), Box<dyn std::error::Error>> {
    let project_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    let python = project_dir.join(".venv").join("bin").join("python");
    let script = project_dir.join("Python").join("train_model.py");

    if !python.exists() {
        return Err(format!("Python virtual environment not found: {}", python.display()).into());
    }

    if !script.exists() {
        return Err(format!("train_model.py not found: {}", script.display()).into());
    }

    println!("Starting ML self-training update...");

    let status = Command::new(&python)
        .arg(&script)
        .arg("update")
        .arg("--epochs")
        .arg("2")
        .current_dir(&project_dir)
        .status()?;

    if status.success() {
        println!("ML model update completed successfully.");
        Ok(())
    } else {
        Err(format!(
            "train_model.py update failed with exit code {:?}",
            status.code()
        )
        .into())
    }
}

pub fn pretrain_ml_model(
    n_examples: usize,
    epochs: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let project_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    let python = project_dir.join(".venv").join("bin").join("python");
    let script = project_dir.join("Python").join("train_model.py");

    eprintln!("Pretrain Python: {}", python.display());
    eprintln!("Python exists: {}", python.exists());
    eprintln!("Pretrain script: {}", script.display());
    eprintln!("Script exists: {}", script.exists());

    if !python.exists() {
        return Err(format!("Python executable does not exist: {}", python.display()).into());
    }

    if !script.exists() {
        return Err(format!("Training script does not exist: {}", script.display()).into());
    }

    let status = Command::new(&python)
        .arg(&script)
        .arg("pretrain")
        .arg("--n")
        .arg(n_examples.to_string())
        .arg("--epochs")
        .arg(epochs.to_string())
        .current_dir(&project_dir)
        .status()
        .map_err(|err| {
            format!(
                "Could not start pretraining.\n\
                 Python: {}\n\
                 Script: {}\n\
                 Underlying error: {}",
                python.display(),
                script.display(),
                err
            )
        })?;

    if !status.success() {
        return Err(format!("Pretraining failed with exit code {:?}", status.code()).into());
    }

    Ok(())
}
