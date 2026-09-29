use std::collections::HashMap;
use std::collections::HashSet;
use std::ffi::{CStr, CString};
use std::os::raw::c_char;
mod ffi;
use ffi::*;

//use librna_sys::*;
use rand::RngExt;
use rand::seq::IndexedRandom;
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{self, Write};
use std::os::raw::c_void;
use std::path::{PathBuf, Path};
use std::process::{Command, Stdio};





const RIBOSOMAL_RNA: bool = false;
// RIBOSOMA_SEQUENCE can be changed to any start sequence of desire. If it is longer than the structure the sequence will be sliced accordingly.
const RIBOSOME_SEQUENCE: &str = "GGCGCCCGCCCCGCGCCGGGGGCCGCCGGCCGGCGGCCCCGGGGGCCCCCCGCCCCCCCGCCCCGCGCGCCCCGCGCGGCGCGGGCCGCGGGCCGGGGGCGGCGCCGGGGCCCCCCCGCCCCCGGCCGGGCCCGGGGCCCGCGCCGGCCCCCCCCGGCGCCGGCGGGCGGGGCCCGCGCCCCGCCCCCGCGCGGCCCGCCCGCCGCGCCCCCCGCCGCCCGGC";
const GC_TEST: bool = false; // This is just if you want your start sequence to be purely paired GC-pairs
const GC_THRESHOLD: f64 = 50.00; // must be a float

fn main() -> io::Result<()> {
    const DEFAULT_N_RUNS: usize = 3;
    const DEFAULT_N_STARTS: i64 = 5;
    const MAX_STEPS: i64 = 2_100;
    const WOBBLE_FREQUENCY: f64 = 0.0;

    println!("========================================");
    println!("RNA design configuration");
    println!("Press Enter to accept a default value.");
    println!("========================================");

    let n_runs = ask_positive_usize(
        "How many complete runs should be performed?",
        DEFAULT_N_RUNS,
    );

    let n_starts = ask_positive_i64(
        "How many hill-climbing starts per run?",
        DEFAULT_N_STARTS,
    );

    println!(
        "\nConfiguration selected: N_RUNS={}, N_STARTS={}, MAX_STEPS={}\n",
        n_runs,
        n_starts,
        MAX_STEPS
    );

    unsafe {
        let mut md: vrna_md_t = std::mem::zeroed();
        vrna_md_set_default(&mut md);
        md.temperature = 37.0;
        md.dangles = 1;

        println!("--- Rust vrna_md_t fields ---");
        println!("min_loop_size={}", md.min_loop_size);
        println!("max_bp_span={}", md.max_bp_span);
        println!("noLP={}", md.noLP);
        println!("noGU={}", md.noGU);
        println!("window_size={}", md.window_size);
    }

    let input_root: PathBuf = if Path::new("misc").is_dir() {
        PathBuf::from("misc")
    } else if Path::new("../misc").is_dir() {
        PathBuf::from("../misc")
    } else {
        panic!(
            "Could not find a misc directory. Looked in:\n\
            - {}\n\
            - {}",
            Path::new("misc").display(),
            Path::new("../misc").display(),
        );
    };

    println!(
        "Current working directory: {}",
        std::env::current_dir()
            .expect("could not determine current directory")
            .display()
    );

    println!("Input directory: {}", input_root.display());

    let input_dir = fs::read_dir(&input_root)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", input_root.display()));

    let mut input_files: Vec<PathBuf> = input_dir
        .filter_map(|entry| {
            let path = entry.ok()?.path();

            let is_input_file = path.is_file()
                && path.file_name()?.to_str()?.starts_with("input_")
                && path.extension()?.to_str()? == "txt";

            is_input_file.then_some(path)
        })
        .collect();

    input_files.sort();

    println!("Found {} input files.", input_files.len());

    if input_files.is_empty() {
        eprintln!(
            "No files matching input_*.txt were found in {}",
            input_root.display()
        );
    }

    let output_base = input_root.join("output");

    fs::create_dir_all(&output_base).unwrap_or_else(|e| {
        panic!(
            "Could not create output directory {}: {e}",
            output_base.display()
        )
    });

    println!("Output base directory: {}", output_base.display());

    
    let mut final_results = String::new();

    for run in 1..=n_runs {
        println!("\n########################################");
        println!("RUN {} of {}", run, n_runs);
        println!("########################################");

        let run_dir = output_base.join(format!("run_{}", run));
        fs::create_dir_all(&run_dir).expect("failed to create run directory");

        for input_path in &input_files {
            let input_path_str = input_path.to_str().expect("Invalid input file path");

            println!("\n========================================");
            println!("PROCESSING: {}", input_path_str);
            println!("========================================");

            let (seq, target_structure) = match read_input_file(input_path_str) {
                Ok((seq, target)) => {
                    println!("Sequence : {}", seq);
                    println!("Structure: {}", target);
                    (seq, target)
                }

                Err(e) => {
                    eprintln!("Failed to read {}: {}", input_path_str, e);
                    continue;
                }
            };

            if RIBOSOMAL_RNA {
                println!("====RIBOSOMAL SEQUENCE USED====");
            }

            if GC_TEST {
                println!("====INITIAL CANDIDATE WILL HAVE OVERLOAD OF GC-PAIRS====");
            }

            let ribo_positions: Vec<usize> = if RIBOSOMAL_RNA {
                seq.char_indices()
                    .filter(|(_, c)| *c == 'N')
                    .map(|(i, _)| i)
                    .collect()
            } else {
                Vec::new()
            };

            let result = if RIBOSOMAL_RNA {
                decomposed_hill_climb_design(
                    &mutate_ks(
                        &insert_ribosome_sequence(&seq),
                        &get_pair_map(&target_structure),
                    ),
                    &target_structure,
                    n_starts,
                    MAX_STEPS,
                    WOBBLE_FREQUENCY,
                    Some(&ribo_positions),
                )
            } else {
                decomposed_hill_climb_design(
                    &mutate_ks(&seq, &get_pair_map(&target_structure)),
                    &target_structure,
                    n_starts,
                    MAX_STEPS,
                    WOBBLE_FREQUENCY,
                    None,
                )
            };

            match result {
                Ok(r) => {
                    
                    println!("\n==== FINAL (run {}) ====", run);
                    println!("sequence     : {}", r.sequence);
                    println!("target       : {}", target_structure);
                    println!("mfe structure: {}", r.mfe_structure);
                    println!("bp_distance  : {}", r.bp_distance);
                    println!("mfe          : {:.2}", r.mfe);
                    println!("Ensemble Diversity: {:.2}", r.ensemble_diversity);
                    println!("slices       : {}", r.n_slices);
                    println!("Ribosomal RNA used: {}", RIBOSOMAL_RNA);
                    

                    match gc_content(&r.sequence) {
                        Some(gc) => println!("GC Content: {:.2}%", gc),
                        None => println!("GC Content: no valid DNA bases found"),
                    }

                    if let Some(identity) = r.ribosome_identity {
                        let n_mismatched = r
                            .ribosome_mismatches
                            .as_ref()
                            .map(|v| v.len())
                            .unwrap_or(0);

                        println!(
                            "ribosome identity  : {:.1}% ({} mismatched positions)",
                            identity * 100.0,
                            n_mismatched
                        );

                        if let Some(mismatches) = &r.ribosome_mismatches {
                            if !mismatches.is_empty() {
                                println!("ribosome mismatches: {:?}", mismatches);
                            }
                        }
                    }

                   
                    final_results.push_str(&format!(
                        "\n========================================\n\
                         INPUT: {}\n\
                         ==== FINAL (run {}) ====\n",
                        input_path_str, run
                    ));

                    final_results.push_str(&format!("sequence     : {}\n", r.sequence));
                    final_results.push_str(&format!("target       : {}\n", target_structure));
                    final_results.push_str(&format!(
                        "mfe structure: {}\n",
                        r.mfe_structure
                    ));
                    final_results.push_str(&format!("bp_distance  : {}\n", r.bp_distance));
                    final_results.push_str(&format!("mfe          : {:.2}\n", r.mfe));
                    final_results.push_str(&format!("Ensemble Diversity: {:.2}\n", r.ensemble_diversity));
                    final_results.push_str(&format!("slices       : {}\n", r.n_slices));
                    final_results.push_str(&format!(
                        "Ribosomal RNA used: {}\n",
                        RIBOSOMAL_RNA
                    ));

                    match gc_content(&r.sequence) {
                        Some(gc) => {
                            final_results.push_str(&format!(
                                "GC Content: {:.2}%\n",
                                gc
                            ));
                        }
                        None => {
                            final_results.push_str(
                                "GC Content: no valid DNA bases found\n"
                            );
                        }
                    }

                    if let Some(identity) = r.ribosome_identity {
                        let n_mismatched = r
                            .ribosome_mismatches
                            .as_ref()
                            .map(|v| v.len())
                            .unwrap_or(0);

                        final_results.push_str(&format!(
                            "ribosome identity  : {:.1}% ({} mismatched positions)\n",
                            identity * 100.0,
                            n_mismatched
                        ));

                        if let Some(mismatches) = &r.ribosome_mismatches {
                            if !mismatches.is_empty() {
                                final_results.push_str(&format!(
                                    "ribosome mismatches: {:?}\n",
                                    mismatches
                                ));
                            }
                        }
                    }

                    
                    let input_filename = input_path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("output");

                    let output_filename = input_filename
                        .strip_prefix("input_")
                        .map(|n| format!("output_{}.txt", n))
                        .unwrap_or_else(|| "output.txt".to_string());

                    let output_path = run_dir.join(&output_filename);

                    let mut file =
                        File::create(&output_path).expect("failed to create output file");

                    writeln!(file, "==== FINAL (run {}) ====", run).unwrap();
                    writeln!(file, "sequence      : {}", r.sequence).unwrap();
                    writeln!(file, "target        : {}", target_structure).unwrap();
                    writeln!(file, "mfe structure : {}", r.mfe_structure).unwrap();
                    writeln!(file, "bp_distance   : {}", r.bp_distance).unwrap();
                    writeln!(file, "mfe           : {:.2}", r.mfe).unwrap();
                    writeln!(file, "Ensemble Diversity : {:.2}", r.ensemble_diversity).unwrap();
                    writeln!(file, "slices        : {}", r.n_slices).unwrap();
                    writeln!(file, "Ribosomal RNA used: {}", RIBOSOMAL_RNA).unwrap();

                    match gc_content(&r.sequence) {
                        Some(gc) => writeln!(file, "GC Content: {:.2}%", gc).unwrap(),
                        None => {
                            writeln!(file, "GC Content: no valid DNA bases found").unwrap()
                        }
                    }

                    if let Some(identity) = r.ribosome_identity {
                        let n_mismatched = r
                            .ribosome_mismatches
                            .as_ref()
                            .map(|v| v.len())
                            .unwrap_or(0);

                        writeln!(
                            file,
                            "ribosome identity  : {:.1}% ({} mismatched positions)",
                            identity * 100.0,
                            n_mismatched
                        )
                        .unwrap();

                        if let Some(mismatches) = &r.ribosome_mismatches {
                            if !mismatches.is_empty() {
                                writeln!(
                                    file,
                                    "ribosome mismatches: {:?}",
                                    mismatches
                                )
                                .unwrap();
                            }
                        }
                    }

                    println!("\nOutput written to {}", output_path.display());
                }

                Err(e) => {
                    eprintln!("Design failed for {} (run {}): {e}", input_path_str, run);

                    // Optional: include failures in the final pager output too.
                    final_results.push_str(&format!(
                        "\n========================================\n\
                         INPUT: {}\n\
                         ==== FAILED (run {}) ====\n\
                         Error: {}\n",
                        input_path_str, run, e
                    ));
                }
            }
        }
    }

    
    if !final_results.is_empty() && ask_to_view_results()? {
        show_in_pager(&final_results)?;
    }

    Ok(())
}

fn get_pair_map(structure: &str) -> HashMap<usize, usize> {
    let mut pair_map = HashMap::new();
    let mut paren_stack: Vec<usize> = Vec::new();
    let mut bracket_stack: Vec<usize> = Vec::new();

    for (i, c) in structure.chars().enumerate() {
        match c {
            '(' => paren_stack.push(i),
            ')' => {
                let j = paren_stack.pop().expect("Unmatched closing parenthesis");
                pair_map.insert(i, j);
                pair_map.insert(j, i);
            }
            '[' => bracket_stack.push(i),
            ']' => {
                let j = bracket_stack.pop().expect("Unmatched closing bracket");
                pair_map.insert(i, j);
                pair_map.insert(j, i);
            }
            _ => {}
        }
    }

    pair_map
}

fn pk_pair_mismatches(seq: &str, pair_map: &HashMap<usize, usize>, target: &str) -> usize {
    let seq_bytes = seq.as_bytes();
    let target_bytes = target.as_bytes();
    let mut mismatches = 0;
    let mut seen = HashSet::new();

    for (&i, &j) in pair_map.iter() {
        if i >= j || seen.contains(&i) {
            continue;
        }

        if target_bytes[i] != b'[' && target_bytes[i] != b']' {
            continue;
        }
        seen.insert(i);
        seen.insert(j);

        let ok = matches!(
            (seq_bytes[i], seq_bytes[j]),
            (b'A', b'U') | (b'U', b'A') | (b'G', b'C') | (b'C', b'G') | (b'G', b'U') | (b'U', b'G')
        );
        if !ok {
            mismatches += 1;
        }
    }
    mismatches
}

fn mutate_seq(seq_in: &str, structure: &str, wobble_frequency: f64, last_global : bool) -> String {
    let mut rng = rand::rng();

    let pair_map = get_pair_map(structure);

    let mut comp_dict = HashMap::new();
    comp_dict.insert('A', 'U');
    comp_dict.insert('U', 'A');
    comp_dict.insert('G', 'C');
    comp_dict.insert('C', 'G');

    let hard_loops: bool = find_hard_loops(structure);

    let _nucleotides = ['A', 'U', 'G', 'C'];

    let paired_nucleotides: &[char] = if GC_TEST || last_global{ // Change here if you want a higher or lower GC-content
            &['G', 'C']
        } else {
            //&['A', 'U', 'G', 'G', 'G', 'C', 'C', 'C', 'G', 'C']
            &['G', 'C'] 
            //&['A', 'U', 'G', 'C']
            //&['A', 'U']
        };

    
    

    

    //let paired_nucleotides = ['A', 'U', 'G', 'G', 'G', 'C', 'C', 'C', 'G', 'C'];

    //let purines = ['A', 'G', 'U'];

    let purines: &[char] = if GC_TEST{
            &['G', 'C']
        } else {
            &['A', 'G', 'U']
        };

    

    let mut mut_seq: Vec<char> = seq_in.chars().collect();
    let mut_struct: Vec<char> = structure.chars().collect();

    for i in 0..mut_seq.len() {
        if mut_seq[i] != 'N' && mut_seq[i] != 'K' && mut_seq[i] != 'S' {
            continue;
        }

        if let Some(&j) = pair_map.get(&i) {
            if mut_seq[j] != 'N' {
                mut_seq[i] = *comp_dict.get(&mut_seq[j]).unwrap();
            } else if mut_seq[i] == 'K' {
                let choice = if rng.random::<f64>() < 0.5 { 'G' } else { 'U' };
                mut_seq[i] = choice;

                if rng.random::<f64>() < wobble_frequency && choice == 'G' {
                    mut_seq[j] = 'U';
                } else {
                    mut_seq[j] = *comp_dict.get(&choice).unwrap();
                }
            } else if mut_seq[i] == 'S' {
                let choice = if rng.random::<f64>() < 0.5 { 'G' } else { 'C' };
                mut_seq[i] = choice;
                mut_seq[j] = *comp_dict.get(&choice).unwrap();
            } else if hard_loops && mut_struct[i] == '(' && mut_struct[i + 1] == ')' {
                let choice = if rng.random::<f64>() < 0.5 { 'G' } else { 'C' };
                mut_seq[i] = choice;
                mut_seq[i + 1] = *comp_dict.get(&choice).unwrap();
                println!("HARD LOOP WARNING -> GC-PAIR!");
                if i + 2 < j {
                    if let Some(&partner_of_next) = pair_map.get(&(i + 2)) {
                        if partner_of_next == j - 2
                            && mut_seq[i + 1] == 'N'
                            && mut_seq[j - 1] == 'N'
                        {
                            let choice2 = if rng.random::<f64>() < 0.5 { 'G' } else { 'C' };
                            mut_seq[i + 2] = choice2;
                            mut_seq[j - 2] = *comp_dict.get(&choice2).unwrap();
                        }
                    }
                }

                if i + 3 < j {
                    if let Some(&partner_of_next) = pair_map.get(&(i + 3)) {
                        if partner_of_next == j - 3
                            && mut_seq[i + 2] == 'N'
                            && mut_seq[j - 2] == 'N'
                        {
                            let choice3 = if rng.random::<f64>() < 0.5 { 'G' } else { 'U' };
                            mut_seq[i + 3] = choice3;
                            mut_seq[j - 3] = *comp_dict.get(&choice3).unwrap();
                        }
                    }
                }
            } else {
                let near_loop =
                    is_adjacent_to_loop(i, structure) || is_adjacent_to_loop(j, structure);

                if near_loop {
                    let choice = if rng.random::<f64>() < 0.5 { 'G' } else { 'C' };
                    mut_seq[i] = choice;
                    mut_seq[j] = *comp_dict.get(&choice).unwrap();

                    if i + 1 < j {
                        if let Some(&partner_of_next) = pair_map.get(&(i + 1)) {
                            if partner_of_next == j - 1
                                && mut_seq[i + 1] == 'N'
                                && mut_seq[j - 1] == 'N'
                            {
                                let choice2 = if rng.random::<f64>() < 0.5 { 'G' } else { 'C' };
                                mut_seq[i + 1] = choice2;
                                mut_seq[j - 1] = *comp_dict.get(&choice2).unwrap();
                            }
                        }
                    }
                    if i + 2 < j {
                        if let Some(&partner_of_next) = pair_map.get(&(i + 2)) {
                            if partner_of_next == j - 2
                                && mut_seq[i + 1] == 'N'
                                && mut_seq[j - 1] == 'N'
                            {
                                let choice2 = if rng.random::<f64>() < 0.5 { 'G' } else { 'C' };
                                mut_seq[i + 2] = choice2;
                                mut_seq[j - 2] = *comp_dict.get(&choice2).unwrap();
                            }
                        }
                    }
                } else {
                    let choice = paired_nucleotides[rng.random_range(0..paired_nucleotides.len())];
                    mut_seq[i] = choice;

                    if rng.random::<f64>() < wobble_frequency && choice == 'G' {
                        mut_seq[j] = 'U';
                    } else {
                        mut_seq[j] = *comp_dict.get(&choice).unwrap();
                    }
                }
            }
        } else {
            mut_seq[i] = purines[rng.random_range(0..purines.len())]; // Hopefully only mutates long unpaired stretches into A's and G's
        }
    }

    mut_seq.into_iter().collect()
}

struct RibosomeSimilarity {
    matches: usize,
    total: usize,
    identity: f64,
    mismatched_positions: Vec<usize>,
}

pub struct GlobalDesignResult {
    sequence: String,
    mfe_structure: String,
    bp_distance: i64,
    mfe: f64,
    n_slices: usize,
    ribosome_identity: Option<f64>,
    ribosome_mismatches: Option<Vec<usize>>,
    ensemble_diversity: f64,         
}
#[derive(Debug, Clone)]
pub struct Substructure {
    pub start: usize,
    pub end: usize,
    pub structure: String,
    pub designable: Vec<usize>,
}

pub struct FoldResult {
    pub bp_distance: i64,
    pub structure: String,
    pub mfe: f64,
}

pub struct DesignResult {
    pub bp_distance: i64,
    pub structure: String,
    pub mfe: f64,
    pub sequence: String,
}

impl Clone for FoldResult {
    fn clone(&self) -> Self {
        FoldResult {
            bp_distance: self.bp_distance,
            structure: self.structure.clone(),
            mfe: self.mfe,
        }
    }
}

impl Clone for DesignResult {
    fn clone(&self) -> Self {
        DesignResult {
            bp_distance: self.bp_distance,
            structure: self.structure.clone(),
            mfe: self.mfe,
            sequence: self.sequence.clone(),
        }
    }
}

pub fn bp_distance_to_target(seq: &str, target: &str) -> DesignResult {
    let n = seq.len();
    let target_no_pk = strip_pseudoknots(target);

    unsafe {
        let mut md: vrna_md_t = std::mem::zeroed();
        vrna_md_set_default(&mut md);
        md.temperature = 37.0;
        md.dangles = 1;

        let seq_bytes = seq.as_bytes();
        let target_bytes = target_no_pk.as_bytes();
        let mut fold_seq_bytes = seq.as_bytes().to_vec();
        for i in 0..n {
            if seq_bytes[i] != b'N' && target_bytes[i] == b'.' {
                fold_seq_bytes[i] = b'N'; // Mask for folding only
            }
        }
        let fold_seq = String::from_utf8(fold_seq_bytes).unwrap();

        /*
        * Pass a borrowed byte slice to CString::new rather than moving fold_seq.
        * This is portable and means fold_seq remains available if needed later.
        */
        let seq_c = CString::new(fold_seq.as_bytes()).expect("seq has interior NUL");

        let fc = vrna_fold_compound(
            seq_c.as_ptr(),
            &md,
            VRNA_OPTION_MFE as u32,
        );

        assert!(!fc.is_null(), "vrna_fold_compound returned null");

        let mut structure: Vec<c_char> = vec![0; n + 1];
        let mfe = vrna_mfe(fc, structure.as_mut_ptr());

        let structure_str = CStr::from_ptr(structure.as_ptr())
            .to_string_lossy()
            .into_owned();

        let target_c = CString::new(target_no_pk.as_bytes())
            .expect("target structure has interior NUL");

        let structure_c = CString::new(structure_str.as_bytes())
            .expect("returned structure has interior NUL");

        let distance = vrna_bp_distance(
            target_c.as_ptr(),
            structure_c.as_ptr(),
        );

        vrna_fold_compound_free(fc);



        DesignResult {
            bp_distance: distance as i64,
            structure: structure_str,
            mfe: mfe as f64,
            sequence: seq.to_string(), // Return the REAL sequence, not the masked one
        }
    }
}

fn hill_climb_design(
    seq_in: &str,
    target_structure: &str,
    n_positions: Vec<usize>,
    max_steps: i64,
    wobble_frequency: f64,
    n_keep: usize,
    initial_temp: f64,
    cooling_rate: f64,
    reheat_after: i64,
    reheat_temp: f64,
    if_slices: bool,
    last_global : bool,
) -> Vec<DesignResult> {
    let pair_map = get_pair_map(target_structure);

    let slice_len = target_structure.len();
    let n_designable = n_positions.len();

    let designable_positions: HashSet<usize> =
        n_positions.iter().copied().collect();
        
    let n_fixed = slice_len.saturating_sub(n_designable);
    let fixed_fraction = n_fixed as f64 / slice_len.max(1) as f64;

    // let dist_threshold: i64 = ((n_fixed as f64) / 12.0).ceil() as i64;

    let dist_threshold: i64 = if RIBOSOMAL_RNA {
        ((n_fixed as f64) / 6.0).ceil() as i64
    } else {
        ((n_fixed as f64) / 12.0).ceil() as i64
    };

    let hard_loops = find_hard_loops(&target_structure);

    let dist_threshold = if hard_loops {
        dist_threshold.max(1) + 2
    } else {
        dist_threshold.max(1)
    };

    let p_threshold: f64 = if fixed_fraction > 0.6 && !last_global {
        0.70
    } else if hard_loops && !last_global{
        0.75
    } else {
        if !last_global{
        0.60 } else {
            0.10
        }
    };

    println!(
        "  [hill_climb] slice_len={} designable={} fixed={} ({:.0}% fixed) → dist_threshold={}, p_threshold={:.2}",
        slice_len,
        n_designable,
        n_fixed,
        fixed_fraction * 100.0,
        dist_threshold,
        p_threshold,
    );

    let nucleotides = vec!['A', 'U', 'G', 'C'];
    let pair_options: Vec<(char, char)> = vec![
        ('A', 'U'),
        ('U', 'A'),
        ('G', 'C'),
        ('C', 'G'),
        ('G', 'U'),
        ('U', 'G'),
    ];
    let seq_in_chars: Vec<char> = seq_in.chars().collect();
    let comp_dict = HashMap::from([('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G')]);

    let mut current = if RIBOSOMAL_RNA {
        mutate_ribosome(seq_in, target_structure, &n_positions)
    } else {
        mutate_seq(seq_in, target_structure, wobble_frequency, last_global)
    };

    //let mut current = mutate_seq(seq_in, target_structure, wobble_frequency);

    let mut current_chars: Vec<char> = current.chars().collect();

    let mut stem_pairs: Vec<(usize, usize)> = n_positions
        .iter()
        .filter_map(|&pos| {
            pair_map.get(&pos).map(|&partner| {
                if pos < partner {
                    (pos, partner)
                } else {
                    (partner, pos)
                }
            })
        })
        .collect();
    stem_pairs.sort_by_key(|&(i, _)| i);
    stem_pairs.dedup();

    for &(pos, partner) in &stem_pairs {
        if seq_in_chars[pos] != 'N' && seq_in_chars[pos] != 'K' && seq_in_chars[pos] != 'S' {
            continue;
        }

        if pos == 0 || partner + 1 >= seq_in_chars.len() {
            continue;
        }
        let outer_5 = pos - 1;
        let outer_3 = partner + 1;

        if pair_map.get(&outer_5) != Some(&outer_3) {
            continue;
        }

        let outer_5_nuc = seq_in_chars[outer_5];
        let outer_3_nuc = seq_in_chars[outer_3];

        if outer_5_nuc == 'N' || outer_3_nuc == 'N' {
            continue;
        }

        let is_weak = matches!(
            (outer_5_nuc, outer_3_nuc),
            ('G', 'U') | ('U', 'G') | ('A', 'U') | ('U', 'A')
        );

        if is_weak {
            current_chars[pos] = 'G';
            current_chars[partner] = 'C';
            println!(
                "  Option B: forcing GC at ({}, {}) — outer pair {}{} is weak",
                pos, partner, outer_5_nuc, outer_3_nuc
            );
        }
    }
    current = current_chars.into_iter().collect();

    let init_result = bp_distance_to_target(&current, target_structure);

    let initial_pk_mismatches =
        pk_pair_mismatches(&current, &pair_map, target_structure);

    let mut current_dist =
        init_result.bp_distance + 2 * initial_pk_mismatches as i64;

    let mut current_structure = init_result.structure.clone();
    let mut current_mfe = init_result.mfe;

    let mut current_energy_gap =
        (energy_of_target_structure(&current, target_structure) - init_result.mfe).max(0.0);

    let mut best_candidates = vec![(
        current_dist,
        current.clone(),
        current_structure.clone(),
        current_mfe,
    )];

    let mut temp = initial_temp;
    let mut rng = rand::rng();
    let mut last_improvement_step: i64 = 0;

    // ################## Temporarily --> This is purely for testing
    let mut n_iterations_for_testing: Vec<i64> = vec![0];
    let mut temperature_for_testing: Vec<f64> = vec![temp];
    let mut current_dist_for_testing: Vec<i64> = vec![current_dist];

    // ##################

    for step in 0..max_steps {
        if best_candidates[0].0 == 0 && best_candidates.len() >= n_keep {
            break;
        }

        if step - last_improvement_step > reheat_after {
            temp = reheat_temp;
            last_improvement_step = step;
            println!("step {step}: REHEATING to temp={temp:.2}");
        }

        let mismatched =
            find_mismatched_positions(&current_structure, target_structure, &n_positions);

        // === THIS IS JUST A TEST CODE. IF THIS NOT WORK GO BACK TO THE OLD WAY =======
        if current_dist > 0 && mismatched.is_empty() && step > 20 {
            let target_no_pk = strip_pseudoknots(target_structure);
            let target_bytes = target_no_pk.as_bytes();
            let current_bytes = current_structure.as_bytes();

            // Indices of *all* mismatches in this slice (local coords)
            let mut all_mismatches = Vec::new();
            for i in 0..current_bytes.len() {
                if current_bytes[i] != target_bytes[i] {
                    all_mismatches.push(i);
                }
            }

            // Fixed mismatched positions = mismatches not in n_positions
            use std::collections::HashSet;
            let designable_set: HashSet<usize> = n_positions.iter().copied().collect();
            let fixed_mismatches: Vec<usize> = all_mismatches
                .into_iter()
                .filter(|i| !designable_set.contains(i))
                .collect();

            if !fixed_mismatches.is_empty() {
                let (_cost, defects, _mfe, _e_target) =
                    compute_pf_defect(&current, target_structure);

                let mut max_p_paired = 0.0;
                for &i in &fixed_mismatches {
                    let p_paired = 1.0 - defects[i];
                    if p_paired > max_p_paired {
                        max_p_paired = p_paired;
                    }
                }

                //let p_threshold = 0.55;
                //let dist_threshold = 3;

                if current_dist <= dist_threshold && max_p_paired < p_threshold {
                    println!(
                        "  early exit: dist={} but remaining mismatches are weakly paired (max p≈{:.2})",
                        current_dist, max_p_paired
                    );
                    break;
                }

                if current_dist <= dist_threshold && if_slices {
                    break;
                };
            }

            println!(
                "  early exit: all mismatches at fixed positions, dist={current_dist} (no PF exception)"
            );
            break;
        }
        // ========= END OF TEST BLOCK ================

        let stuck = step - last_improvement_step > 100;
        let use_double = stuck && current_dist <= 2;

        let mut candidate: Vec<char> = current.chars().collect();

        let target_bytes = target_structure.as_bytes();
        let stuck_positions: Vec<usize> = n_positions
            .iter()
            .filter(|&&p| {
                p < current_structure.len() && current_structure.as_bytes()[p] != target_bytes[p]
            })
            .copied()
            .collect();

                if use_double {
            let mut search_positions: Vec<usize> = if stuck_positions.len() >= 2 {
                stuck_positions.clone()
            } else {
                n_positions.clone()
            };

            if search_positions.len() > 30 {
                use rand::seq::SliceRandom;

                search_positions.shuffle(&mut rng);
                search_positions.truncate(30);
            }

            let mut best_trial = candidate.clone();
            let mut best_cost = f64::MAX;

            for idx1 in 0..search_positions.len() {
                for idx2 in (idx1 + 1)..search_positions.len() {
                    let i1 = search_positions[idx1];
                    let i2 = search_positions[idx2];

                    
                    if pair_map.get(&i1) == Some(&i2) {
                        continue;
                    }

                    let opts1: Vec<char> = if let Some(&j1) = pair_map.get(&i1) {
                        if designable_positions.contains(&j1) {
                            nucleotides.clone()
                        } else {
                            let fixed_partner = candidate[j1];

                            vec![
                                *comp_dict
                                    .get(&fixed_partner)
                                    .expect(
                                        "fixed paired position must contain A, U, G, or C",
                                    ),
                            ]
                        }
                    } else {
                        nucleotides.clone()
                    };

                    let opts2: Vec<char> = if let Some(&j2) = pair_map.get(&i2) {
                        if designable_positions.contains(&j2) {
                            nucleotides.clone()
                        } else {
                            let fixed_partner = candidate[j2];

                            vec![
                                *comp_dict
                                    .get(&fixed_partner)
                                    .expect(
                                        "fixed paired position must contain A, U, G, or C",
                                    ),
                            ]
                        }
                    } else {
                        nucleotides.clone()
                    };

                    for &n1 in &opts1 {
                        for &n2 in &opts2 {
                            let mut trial = candidate.clone();

                            
                            trial[i1] = n1;

                            if let Some(&j1) = pair_map.get(&i1) {
                                if designable_positions.contains(&j1) {
                                    trial[j1] = comp_dict[&n1];
                                }
                            }

                            
                            trial[i2] = n2;

                            if let Some(&j2) = pair_map.get(&i2) {
                                if designable_positions.contains(&j2) {
                                    trial[j2] = comp_dict[&n2];
                                }
                            }

                            let trial_str: String = trial.iter().collect();

                            let trial_result =
                                bp_distance_to_target(&trial_str, target_structure);

                            let pk_mismatches =
                                pk_pair_mismatches(&trial_str, &pair_map, target_structure);

                            let trial_cost =
                                trial_result.bp_distance as f64
                                    + 2.0 * pk_mismatches as f64;

                            if trial_cost < best_cost {
                                best_cost = trial_cost;
                                best_trial = trial;
                            }
                        }
                    }
                }
            }

            candidate = best_trial;
        } else {
            // === SINGLE MUTATION ===
            let i: usize = if !mismatched.is_empty() && rng.random::<f64>() < 0.7 {
                *mismatched.choose(&mut rng).unwrap()
            } else {
                *n_positions.choose(&mut rng).unwrap()
            };

            if let Some(&j) = pair_map.get(&i) {
                
                if !designable_positions.contains(&j) {
                    let current_nuc = candidate[i];
                    let mut best_nuc = current_nuc;
                    let mut best_cost = f64::MAX;

                    
                    let fixed_partner = candidate[j];

                    let allowed_nucs = vec![
                        *comp_dict
                            .get(&fixed_partner)
                            .expect("fixed paired position must contain A, U, G, or C"),
                    ];

                    for &nuc in &nucleotides {
                        let mut trial = candidate.clone();
                        trial[i] = nuc;

                        let trial_str: String = trial.iter().collect();

                        let trial_result =
                            bp_distance_to_target(&trial_str, target_structure);

                        let pk_mismatches =
                            pk_pair_mismatches(&trial_str, &pair_map, target_structure);

                        let trial_cost =
                            trial_result.bp_distance as f64
                                + 2.0 * pk_mismatches as f64;

                        if trial_cost < best_cost {
                            best_cost = trial_cost;
                            best_nuc = nuc;
                        }
                    }

                    if best_nuc == current_nuc {
                        best_nuc = allowed_nucs[0];
                    }

                    candidate[i] = best_nuc;
                } else {
                    
                    let current_pair = (candidate[i], candidate[j]);
                    let mut best_pair = current_pair;
                    let mut best_cost = f64::MAX;

                    for &(ni, nj) in &pair_options {
                        let mut trial = candidate.clone();
                        trial[i] = ni;
                        trial[j] = nj;

                        let trial_str: String = trial.iter().collect();

                        let trial_result =
                            bp_distance_to_target(&trial_str, target_structure);

                        let pk_mismatches =
                            pk_pair_mismatches(&trial_str, &pair_map, target_structure);

                        let trial_cost =
                            trial_result.bp_distance as f64
                                + 2.0 * pk_mismatches as f64;

                        if trial_cost < best_cost {
                            best_cost = trial_cost;
                            best_pair = (ni, nj);
                        }
                    }

                    if best_pair == current_pair {
                        best_pair =
                            pair_options[rng.random_range(0..pair_options.len())];
                    }

                    candidate[i] = best_pair.0;
                    candidate[j] = best_pair.1;
                }
            } else {
                
                let current_nuc = candidate[i];
                let mut best_nuc = current_nuc;
                let mut best_cost = f64::MAX;

                for &nuc in &nucleotides {
                    let mut trial = candidate.clone();
                    trial[i] = nuc;

                    let trial_str: String = trial.iter().collect();

                    let trial_result =
                        bp_distance_to_target(&trial_str, target_structure);

                    let pk_mismatches =
                        pk_pair_mismatches(&trial_str, &pair_map, target_structure);

                    let trial_cost =
                        trial_result.bp_distance as f64
                            + 2.0 * pk_mismatches as f64;

                    if trial_cost < best_cost {
                        best_cost = trial_cost;
                        best_nuc = nuc;
                    }
                }

                if best_nuc == current_nuc {
                    best_nuc = nucleotides[rng.random_range(0..nucleotides.len())];
                }

                candidate[i] = best_nuc;
            }
        }

        let candidate_string: String = candidate.into_iter().collect();


        
        let guard_result = bp_distance_to_target(&candidate_string, target_structure);
        let cand_pk_mismatches = pk_pair_mismatches(&candidate_string, &pair_map, target_structure);
        let cand_dist = guard_result.bp_distance + 2 * cand_pk_mismatches as i64;
        let cand_energy_gap = (energy_of_target_structure(&candidate_string, target_structure)
            - guard_result.mfe)
            .max(0.0);

        let accept = if cand_dist < current_dist {
            true
        } else if cand_dist == current_dist {
            cand_energy_gap < current_energy_gap
                || rng.random::<f64>()
                    < (-(cand_energy_gap - current_energy_gap) / temp.max(1e-6)).exp()
        } else {
            rng.random::<f64>() < (-(cand_dist - current_dist) as f64 / temp.max(1e-6)).exp()
        };

        n_iterations_for_testing.push(step);
        temperature_for_testing.push(temp);
        current_dist_for_testing.push(current_dist);

        let current_bytes = current_structure.as_bytes();
        let mut all_mismatches = Vec::new();
        for i in 0..current_bytes.len() {
            if current_bytes[i] != target_bytes[i] {
                all_mismatches.push(i);
            }
        }

        use std::collections::HashSet;
        let designable_set: HashSet<usize> = n_positions.iter().copied().collect();
        let fixed_mismatches: Vec<usize> = all_mismatches
            .into_iter()
            .filter(|i| !designable_set.contains(i))
            .collect();

        if !fixed_mismatches.is_empty() {
            let (_cost, defects, _mfe, _e_target) = compute_pf_defect(&current, target_structure);

            let mut max_p_paired = 0.0;
            for &i in &fixed_mismatches {
                let p_paired = 1.0 - defects[i];
                if p_paired > max_p_paired {
                    max_p_paired = p_paired;
                }
            }

            //let p_threshold = 0.55;
            //let dist_threshold = 3;

            if cand_dist <= dist_threshold && max_p_paired < p_threshold {
                println!(
                    "  early exit: dist={} but remaining mismatches are weakly paired (max p≈{:.2})",
                    cand_dist, max_p_paired
                );
                break;
            }
        }

        if cand_dist <= dist_threshold && if_slices {
            println!(
                "Early exit for slices at dist={}: Remaining positions will be solved globally",
                cand_dist
            );
            break;
        };

        // End of break block

        if accept {
            let improved = cand_dist < current_dist;
            current = candidate_string.clone();
            current_dist = cand_dist;

            current_energy_gap = (energy_of_target_structure(&candidate_string, target_structure)
                - guard_result.mfe)
                .max(0.0);

            current_structure = guard_result.structure.clone();
            current_mfe = guard_result.mfe;

            if improved {
                last_improvement_step = step;
            }

            let existing_seqs: HashSet<&String> = best_candidates.iter().map(|c| &c.1).collect();
            if !existing_seqs.contains(&candidate_string) {
                best_candidates.push((
                    current_dist,
                    candidate_string,
                    current_structure.clone(),
                    current_mfe,
                ));
                best_candidates.sort_by(|a, b| a.0.cmp(&b.0));
                best_candidates.truncate(n_keep);
            }
        }

        temp *= cooling_rate;
        temp = temp.max(0.05);

        if step % 500 == 0 {
            println!(
                "step {step}: dist={current_dist}, gap={current_energy_gap:.2}, best_dist={}, temp={temp:.4}, mode={}",
                best_candidates[0].0,
                if use_double { "DOUBLE" } else { "single" }
            );
        }
    }

    best_candidates
        .into_iter()
        .map(|(bp_distance, sequence, structure, mfe)| DesignResult {
            bp_distance,
            structure,
            mfe,
            sequence,
        })
        .collect()
}
    

fn multi_start_hill_climb_design(
    seq_in: &str,
    target_structure: &str,
    n_position: Vec<usize>,
    n_runs: i64,
    max_steps: i64,
    wobble_frequency: f64,
    if_slices: bool,
    last_global : bool,
) -> Vec<DesignResult> {
    let mut all_candidates: Vec<DesignResult> = Vec::new();
    let n_keep_per_run: usize = 1;
    let temp: f64 = 1.0;
    let final_temp: f64 = 0.01;
    let cooling_rate: f64 = (final_temp / temp).powf(1.0 / max_steps as f64);
    //let cooling_rate: f64 = 0.999; // <-- optional override for testing
    let reheat_after: i64 = 500; // number of stagnant steps before reheating
    let reheat_temp: f64 = 1.0; // temperature to jump back

    let all_pools: Vec<Vec<DesignResult>> = (0..n_runs)
        .into_par_iter()
        .map(|_run| {
            hill_climb_design(
                seq_in,
                target_structure,
                n_position.clone(),
                max_steps,
                wobble_frequency,
                n_keep_per_run,
                temp,
                cooling_rate,
                reheat_after,
                reheat_temp,
                if_slices,
                last_global,
            )
        })
        .collect();

    for pool in all_pools {
        all_candidates.extend(pool);
    }

    let _n_keep_global: usize = 20;
    all_candidates.sort_by(|a, b| {
        a.bp_distance.cmp(&b.bp_distance).then_with(|| {
            a.mfe
                .partial_cmp(&b.mfe)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    });
    all_candidates
}

fn get_pk_pairs(structure: &str) -> Vec<(usize, usize)> {
    let pair_map = get_pair_map(structure);
    let b = structure.as_bytes();
    let mut pk_pairs = Vec::new();
    let mut seen = HashSet::new();
    for (&i, &j) in pair_map.iter() {
        if i >= j || seen.contains(&i) {
            continue;
        }
        if b[i] == b'[' || b[i] == b']' {
            seen.insert(i);
            seen.insert(j);
            pk_pairs.push((i.min(j), i.max(j)));
        }
    }
    pk_pairs
}

pub fn decomposed_hill_climb_design(
    start_seq: &str,
    target: &str,
    n_starts: i64,
    max_steps: i64,
    wobble_frequency: f64,
    ribo_positions: Option<&[usize]>,
) -> Result<GlobalDesignResult, String> {
    assert_eq!(
        start_seq.len(),
        target.len(),
        "start_seq and target must have equal length"
    );

    let ribo_set: HashSet<usize> = ribo_positions
        .map(|p| p.iter().copied().collect())
        .unwrap_or_default();
    if RIBOSOMAL_RNA {
        println!(
            "Restricting mutations to {} ribosomal-inserted positions",
            ribo_set.len()
        );
    }

    let conserved_motifs = important_motifs();
    let conserved_motifs = find_conserved_positions(&start_seq, &conserved_motifs);

    let if_slices: bool = true;
    let slices = decompose(target).map_err(|e| e.to_string())?;

    let mut seq = start_seq.as_bytes().to_vec();
    if RIBOSOMAL_RNA {
        let mut rng = rand::rng();
        let nucleotides = ['A', 'U', 'G', 'C'];
        for pos in 0..seq.len() {
            if !ribo_set.contains(&pos) && matches!(seq[pos], b'N' | b'K' | b'S') {
                seq[pos] = match seq[pos] {
                    b'K' => {
                        if rng.random::<f64>() < 0.5 {
                            b'G'
                        } else {
                            b'U'
                        }
                    }
                    b'S' => {
                        if rng.random::<f64>() < 0.5 {
                            b'G'
                        } else {
                            b'C'
                        }
                    }
                    _ => nucleotides[rng.random_range(0..4)] as u8,
                };
            }
        }
    }

    println!(
        "decomposed {} slices (post-order: children before parents)",
        slices.len()
    );

    for (idx, sub) in slices.iter().enumerate() {
        let local_len = sub.end - sub.start;
        let local_protected_positions =
            get_local_protected_positions(&conserved_motifs, sub.start, sub.end);
        let local_designable: Vec<usize> = if RIBOSOMAL_RNA {
            sub.designable
                .iter()
                .filter(|&&g| ribo_set.contains(&g))
                .map(|&g| g - sub.start)
                .filter(|local_pos| !local_protected_positions.contains(local_pos))
                .collect()
        } else {
            sub.designable
                .iter()
                .filter(|&&g| matches!(start_seq.as_bytes()[g], b'N' | b'K' | b'S'))
                .map(|&g| g - sub.start)
                .collect()
        };

        if local_designable.is_empty() {
            println!(
                "[slice {}/{}] range=[{}, {}] — no designable positions, skipping",
                idx + 1,
                slices.len(),
                sub.start,
                sub.end
            );
            continue;
        }

        let mut slice_start_bytes = seq[sub.start..sub.end].to_vec();
        for &loc in &local_designable {
            slice_start_bytes[loc] = b'N';
        }
        let slice_start = String::from_utf8(slice_start_bytes).unwrap();

        println!(
            "[slice {}/{}] range=[{}, {}] len={} designable={:?} (local={:?})",
            idx + 1,
            slices.len(),
            sub.start,
            sub.end,
            local_len,
            sub.designable,
            local_designable,
        );
        println!("  target:  {}", sub.structure);
        println!("  input:   {}", slice_start);

        let n_designable = local_designable.len();
        let n_total = sub.end - sub.start;

        let pool = multi_start_hill_climb_design(
            &slice_start,
            &sub.structure,
            local_designable,
            n_starts,
            max_steps,
            wobble_frequency,
            if_slices,
            false,
        );

        let best = pool
            .into_iter()
            .min_by_key(|r| r.bp_distance)
            .ok_or_else(|| format!("slice {} returned empty pool", idx))?;

        let fixed_fraction = 1.0 - (n_designable as f64 / n_total as f64);

        if best.bp_distance > 0 && fixed_fraction > 0.5 {
            println!(
                "  NOTE: slice {} best dist={} ({}% fixed), proceeding to next slice",
                idx + 1,
                best.bp_distance,
                (fixed_fraction * 100.0) as u32
            );
        }

        println!(
            "  best:    bp_distance={}, mfe={:.2}",
            best.bp_distance, best.mfe,
        );
        println!("  seq:     {}", best.sequence);
        println!("  struct:  {}", best.structure);

        let designed_bytes = best.sequence.as_bytes();

        if designed_bytes.len() != local_len {
            return Err(format!(
                "slice {} returned sequence of length {}, expected {}",
                idx,
                designed_bytes.len(),
                local_len,
            ));
        }

        seq[sub.start..sub.end].copy_from_slice(designed_bytes);
    }

    let mut rng = rand::rng();
    let nucleotides = ['A', 'U', 'G', 'C'];
    let mut full_seq_chars: Vec<char> = String::from_utf8(seq.clone()).unwrap().chars().collect();

    for c in full_seq_chars.iter_mut() {
        if *c == 'N' {
            *c = nucleotides[rng.random_range(0..nucleotides.len())];
        }
        if *c == 'K' {
            *c = if rng.random::<f64>() < 0.5 { 'G' } else { 'U' };
        }
        if *c == 'S' {
            *c = if rng.random::<f64>() < 0.5 { 'G' } else { 'C' };
        }
    }

    let pk_pairs = get_pk_pairs(target);
    if !pk_pairs.is_empty() {
        let valid_pairs = [
            ('A', 'U'),
            ('U', 'A'),
            ('G', 'C'),
            ('C', 'G'),
            ('G', 'U'),
            ('U', 'G'),
        ];
        let mut n_pk_fixed = 0;
        for &(i, j) in &pk_pairs {
            if RIBOSOMAL_RNA && !ribo_set.contains(&i) && !ribo_set.contains(&j) {
                continue;
            }

            let (bi, bj) = (full_seq_chars[i], full_seq_chars[j]);
            if !valid_pairs.contains(&(bi, bj)) {
                let (ni, nj) = valid_pairs[rng.random_range(0..valid_pairs.len())];
                full_seq_chars[i] = ni;
                full_seq_chars[j] = nj;
                n_pk_fixed += 1;
            }
        }
        if n_pk_fixed > 0 {
            println!("PK repair: fixed {} pseudoknotted pair(s)", n_pk_fixed);
        }
    }

    let full_seq = full_seq_chars.into_iter().collect::<String>();

    println!("Proceeding to global check");
    let pre_repair_result = bp_distance_to_target(&full_seq, target);

    let full_seq = if pre_repair_result.bp_distance > 0 {
        let mismatched_positions: Vec<usize> = if RIBOSOMAL_RNA {
            identify_mismatches(&pre_repair_result.structure, target)
                .into_iter()
                .filter(|&g| ribo_set.contains(&g))
                .collect()
        } else {
            identify_mismatches(&pre_repair_result.structure, target)
        };
        println!(
            "attempting global repair via multi_start_hill_climb_design, bp_distance={}, {} mismatched positions",
            pre_repair_result.bp_distance,
            mismatched_positions.len()
        );

        let repair_max_steps: i64 =
            if ((pre_repair_result.bp_distance as f64 / 0.008).round() as i64) < 800 {
                (pre_repair_result.bp_distance as f64 / 0.008).round() as i64
            } else {
                800
            };
        println!(
            "Max number of iterations for global repair: {}",
            repair_max_steps
        );

        let if_slices: bool = false;

        let repair_pool = multi_start_hill_climb_design(
            &full_seq,
            target,
            mismatched_positions,
            n_starts,
            repair_max_steps,
            wobble_frequency,
            if_slices,
            false,
        );

        match repair_pool.into_iter().min_by_key(|r| r.bp_distance) {
            Some(repaired) if repaired.bp_distance < pre_repair_result.bp_distance => {
                println!(
                    "global repair improved: {} -> {}",
                    pre_repair_result.bp_distance, repaired.bp_distance
                );
                repaired.sequence
            }
            _ => {
                println!("global repair made no improvement, keeping pre-repair sequence");
                full_seq
            }
        }
    } else {
        full_seq
    };

    
    let final_check_result = bp_distance_to_target(&full_seq, target);

    let full_seq = if final_check_result.bp_distance > 0 {
        let remaining_mismatches: Vec<usize> = if RIBOSOMAL_RNA {
            identify_mismatches(&final_check_result.structure, target)
                .into_iter()
                .filter(|&g| ribo_set.contains(&g))
                .collect()
        } else {
            identify_mismatches(&final_check_result.structure, target)
        };

        if remaining_mismatches.is_empty() {
            full_seq
        } else {
            
            let final_repair_max_steps: i64 =
                ((final_check_result.bp_distance as f64 / 0.008).round() as i64)
                    .clamp(1, 100);

            println!(
                "attempting final focused global repair via multi_start_hill_climb_design, \
                 bp_distance={}, {} remaining mismatched positions, max_steps={}",
                final_check_result.bp_distance,
                remaining_mismatches.len(),
                final_repair_max_steps
            );

            let if_slices = false;

            let final_repair_pool = multi_start_hill_climb_design(
                &full_seq,
                target,
                remaining_mismatches,
                n_starts,
                final_repair_max_steps,
                wobble_frequency,
                if_slices,
                true,
            );

            match final_repair_pool.into_iter().min_by_key(|r| r.bp_distance) {
                Some(repaired) if repaired.bp_distance < final_check_result.bp_distance => {
                    println!(
                        "final focused global repair improved: {} -> {}",
                        final_check_result.bp_distance,
                        repaired.bp_distance
                    );
                    repaired.sequence
                }
                _ => {
                    println!(
                        "final focused global repair made no improvement, keeping current sequence"
                    );
                    full_seq
                }
            }
        }
    } else {
        full_seq
    };

    let final_check_result = bp_distance_to_target(&full_seq, target);

    let full_seq = if final_check_result.bp_distance > 0 {
        let remaining_mismatches: Vec<usize> = if RIBOSOMAL_RNA {
            identify_mismatches(&final_check_result.structure, target)
                .into_iter()
                .filter(|&g| ribo_set.contains(&g))
                .collect()
        } else {
            identify_mismatches(&final_check_result.structure, target)
        };

        if remaining_mismatches.is_empty() {
            full_seq
        } else {
            
            let final_repair_max_steps: i64 =
                ((final_check_result.bp_distance as f64 / 0.008).round() as i64)
                    .clamp(1, 100);

            println!(
                "attempting a second final focused global repair via multi_start_hill_climb_design, \
                 bp_distance={}, {} remaining mismatched positions, max_steps={}",
                final_check_result.bp_distance,
                remaining_mismatches.len(),
                final_repair_max_steps
            );

            let if_slices = false;

            let final_repair_pool = multi_start_hill_climb_design(
                &full_seq,
                target,
                remaining_mismatches,
                n_starts,
                final_repair_max_steps,
                wobble_frequency,
                if_slices,
                true,
            );

            match final_repair_pool.into_iter().min_by_key(|r| r.bp_distance) {
                Some(repaired) if repaired.bp_distance < final_check_result.bp_distance => {
                    println!(
                        "the second final focused global repair improved: {} -> {}",
                        final_check_result.bp_distance,
                        repaired.bp_distance
                    );
                    repaired.sequence
                }
                _ => {
                    println!(
                        "the second final focused global repair made no improvement, keeping current sequence"
                    );
                    full_seq
                }
            }
        }
    } else {
        full_seq
    };

    let full_seq = if RIBOSOMAL_RNA || GC_TEST {
        println!("Skipping GC cleanup");
        full_seq
    } else {
        gc_cleanup(&full_seq, start_seq, target, GC_THRESHOLD)
    };

    let result = bp_distance_to_target(&full_seq, target);

    let (mfe_struct_with_pk, n_pk_reinserted) = annotate_pk(&full_seq, target, &result.structure);

    let result_no_pk = bp_distance_to_target(&full_seq, target);
    let pair_map = get_pair_map(target);
    let pk_penalty = pk_pair_mismatches(&full_seq, &pair_map, target) as i64;

    let ensemble_diversity_score: f64 = ensemble_diversity(&full_seq, target);

    println!("---- global verification ----");
    println!("designed:                 {}", full_seq);
    println!("target:                   {}", target);
    println!("mfe struct (VRNA-only):   {}", result.structure);
    println!("mfe struct (with pk):     {}", mfe_struct_with_pk);
    println!("pk pairs reinserted:      {}", n_pk_reinserted);
    println!("bp_distance (VRNA only):  {}", result_no_pk.bp_distance);
    println!("pk penalty (bad pk pairs): {}", pk_penalty);
    println!("bp_distance (pk-aware):   {}", result.bp_distance);
    println!("mfe (VRNA, no pk energy): {:.2}", result.mfe);
    println!("ensemble diversity:       {:.2}", ensemble_diversity_score);

    let pk_predictions = fold_with_pkplex(&full_seq);
    if !pk_predictions.is_empty() {
        println!("---- pkplex verification ----");
        for (k, pk) in pk_predictions.iter().enumerate() {
            println!(
                "  pk {}: dGpk={:.2}, dG1={:.2}, dG2={:.2}, dGint={:.2}, range 5'=[{},{}], 3'=[{},{}]",
                k, pk.dgpk, pk.dg1, pk.dg2, pk.dgint, pk.start_5, pk.end_5, pk.start_3, pk.end_3
            );
        }
    } else {
        println!("---- pkplex: no pseudoknot predicted ----");
    }

    let (ribosome_identity, ribosome_mismatches) = if RIBOSOMAL_RNA {
        let ribo_positions_vec: Vec<usize> = ribo_positions.unwrap_or(&[]).to_vec();
        let sim = ribosome_similarity(&full_seq, &ribo_positions_vec);
        println!(
            "ribosome identity: {}/{} ({:.1}%)",
            sim.matches,
            sim.total,
            sim.identity * 100.0
        );
        (Some(sim.identity), Some(sim.mismatched_positions))
    } else {
        (None, None)
    };

    Ok(GlobalDesignResult {
        sequence: full_seq,
        mfe_structure: mfe_struct_with_pk,
        bp_distance: result.bp_distance,
        mfe: result.mfe,
        n_slices: slices.len(),
        ribosome_identity,
        ribosome_mismatches,
        ensemble_diversity: ensemble_diversity_score,   
    })
}

fn compute_partners(structure: &str) -> Option<Vec<Option<usize>>> {
    let b = structure.as_bytes();
    let n = b.len();
    let mut partner = vec![None; n];
    let mut paren_stack = Vec::new();
    let mut bracket_stack = Vec::new();

    for (i, &c) in b.iter().enumerate() {
        match c {
            b'(' => paren_stack.push(i),
            b')' => {
                let j = paren_stack.pop()?;
                partner[i] = Some(j);
                partner[j] = Some(i);
            }
            b'[' => bracket_stack.push(i),
            b']' => {
                let j = bracket_stack.pop()?;
                partner[i] = Some(j);
                partner[j] = Some(i);
            }
            b'.' | b'_' => {}
            _ => return None,
        }
    }
    if paren_stack.is_empty() && bracket_stack.is_empty() {
        Some(partner)
    } else {
        None
    }
}

fn find_blocks(partner: &[Option<usize>], b: &[u8]) -> Vec<(usize, usize)> {
    let n = partner.len();
    let mut blocks = Vec::new();
    let mut i = 0;
    while i < n {
        let mut reach = i;
        let mut k = i;
        loop {
            if b[k] != b'[' && b[k] != b']' {
                if let Some(p) = partner[k] {
                    if p > reach {
                        reach = p;
                    }
                }
            }
            if k == reach {
                break;
            }
            k += 1;
        }
        blocks.push((i, reach));
        i = reach + 1;
    }
    blocks
}

pub fn decompose(structure: &str) -> Result<Vec<Substructure>, &'static str> {
    let partner = compute_partners(structure).ok_or("unbalanced or unsupported structure")?;
    let b = structure.as_bytes();
    let mut out = Vec::new();

    for (start, reach) in find_blocks(&partner, b) {
        if start == reach {
            continue; // single unpaired position
        }
        if b[start] == b'(' && partner[start] == Some(reach) {
            process_stem_loop(b, &partner, start, reach, &mut out);
        } else {
            let designable: Vec<usize> = (start..=reach).collect();
            out.push(Substructure {
                start,
                end: reach + 1,
                structure: std::str::from_utf8(&b[start..=reach]).unwrap().to_string(),
                designable,
            });
        }
    }

    for sub in out.iter_mut() {
        sub.structure = strip_pseudoknots(&sub.structure);
    }

    assert_slices_balanced(&out).map_err(|e| {
        eprintln!("{e}");
        "unbalanced slice produced by decompose"
    })?;

    Ok(out)
}

pub fn decompose_4(structure: &str) -> Result<Vec<Substructure>, &'static str> {
    let partner = compute_partners(structure).ok_or("unbalanced or unsupported structure")?;
    let b = structure.as_bytes();
    let mut out = Vec::new();

    for (start, reach) in find_blocks(&partner, b) {
        if start == reach {
            continue; // single unpaired position, nothing to do
        }
        if (b[start] == b'(' || b[start] == b'[') && partner[start] == Some(reach) {
            process_stem_loop(b, &partner, start, reach, &mut out);
        } else {
            let designable: Vec<usize> = (start..=reach).collect();
            out.push(Substructure {
                start,
                end: reach + 1,
                structure: std::str::from_utf8(&b[start..=reach]).unwrap().to_string(),
                designable,
            });
        }
    }

    Ok(out)
}

pub fn decompose3(structure: &str) -> Result<Vec<Substructure>, &'static str> {
    let partner = compute_partners(structure).ok_or("unbalanced or unsupported structure")?;
    let b = structure.as_bytes();
    let n = b.len();
    let mut out = Vec::new();
    let mut i = 0;
    while i < n {
        if b[i] == b'(' {
            let j = partner[i].unwrap();
            process_stem_loop(b, &partner, i, j, &mut out);
            i = j + 1;
        } else {
            i += 1;
        }
    }
    Ok(out)
}

pub fn decompose2(structure: &str) -> Result<Vec<Substructure>, &'static str> {
    let n = structure.len();

    compute_partners(structure).ok_or("unbalanced or unsupported structure")?;

    Ok(vec![Substructure {
        start: 0,
        end: n,
        structure: structure.to_string(),
        designable: (0..n).collect(),
    }])
}

fn _process_stem_loop(
    b: &[u8],
    partner: &[Option<usize>],
    i: usize,
    j: usize,
    out: &mut Vec<Substructure>,
) {
    let mut i_inner = i;
    let mut j_inner = j;
    while i_inner + 1 < j_inner
        && b[i_inner + 1] == b'('
        && b[j_inner - 1] == b')'
        && partner[i_inner + 1] == Some(j_inner - 1)
    {
        i_inner += 1;
        j_inner -= 1;
    }

    let mut designable = Vec::new();
    for k in i..=i_inner {
        designable.push(k);
    }
    for k in j_inner..=j {
        designable.push(k);
    }

    let mut k = i_inner + 1;
    while k < j_inner {
        match b[k] {
            b'(' => {
                let l = partner[k].expect("unmatched '(' inside slice");
                process_stem_loop(b, partner, k, l, out);
                k = l + 1;
            }
            b'[' => {
                if let Some(l) = partner[k] {
                    let lo = k.min(l);
                    let hi = k.max(l);
                    let designable_pk: Vec<usize> = (lo..=hi).collect();
                    out.push(Substructure {
                        start: lo,
                        end: hi + 1,
                        structure: std::str::from_utf8(&b[lo..=hi]).unwrap().to_string(),
                        designable: designable_pk,
                    });
                    k = hi + 1;
                } else {
                    designable.push(k);
                    k += 1;
                }
            }

            b')' => {
                designable.push(k);
                k += 1;
            }
            _ => {
                designable.push(k);
                k += 1;
            }
        }
    }

    out.push(Substructure {
        start: i,
        end: j + 1,
        structure: std::str::from_utf8(&b[i..=j]).unwrap().to_string(),
        designable,
    });
}

fn process_stem_loop(
    b: &[u8],
    partner: &[Option<usize>],
    i: usize,
    j: usize,
    out: &mut Vec<Substructure>,
) {
    let mut i_inner = i;
    let mut j_inner = j;
    while i_inner + 1 < j_inner && partner[i_inner + 1] == Some(j_inner - 1) {
        i_inner += 1;
        j_inner -= 1;
    }

    let mut designable = Vec::new();

    for k in i..=i_inner {
        designable.push(k);
    }

    for k in j_inner..=j {
        designable.push(k);
    }

    let mut k = i_inner + 1;
    while k < j_inner {
        if b[k] == b'(' {
            let l = partner[k].unwrap();
            process_stem_loop(b, partner, k, l, out);
            k = l + 1;
        } else {
            designable.push(k);
            k += 1;
        }
    }

    out.push(Substructure {
        start: i,
        end: j + 1,
        structure: std::str::from_utf8(&b[i..=j]).unwrap().to_string(),
        designable,
    });
}

fn energy_of_target_structure(seq: &str, target: &str) -> f64 {
    let target_no_pk = strip_pseudoknots(target);
    unsafe {
        let mut md: vrna_md_t = std::mem::zeroed();
        vrna_md_set_default(&mut md);
        md.temperature = 37.0;
        md.dangles = 1;

        let seq_c = CString::new(seq).unwrap();
        let fc = vrna_fold_compound(seq_c.as_ptr(), &md, VRNA_OPTION_MFE as u32);
        let target_c = CString::new(target_no_pk).unwrap();
        let energy = vrna_eval_structure(fc, target_c.as_ptr());
        vrna_fold_compound_free(fc);
        energy as f64
    }
}

fn find_mismatched_positions(
    current_structure: &str,
    target: &str,
    n_positions: &[usize],
) -> Vec<usize> {
    let current_bytes = current_structure.as_bytes();
    let tgt_no_pk = strip_pseudoknots(target);
    let target_bytes = tgt_no_pk.as_bytes();
    n_positions
        .iter()
        .filter(|&&i| target_bytes[i] != b'.' && current_bytes[i] != target_bytes[i])
        .copied()
        .collect()
}

fn compute_pf_defect(seq: &str, target: &str) -> (f64, Vec<f64>, f64, f64) {
    let n = seq.len();
    let target_bytes = target.as_bytes();
    let target_no_pk = strip_pseudoknots(target);

    unsafe {
        let mut md: vrna_md_t = std::mem::zeroed();
        vrna_md_set_default(&mut md);
        md.temperature = 37.0;
        md.dangles = 1;

        let seq_c = CString::new(seq).expect("seq has interior NUL");
        let fc = vrna_fold_compound(
            seq_c.as_ptr(),
            &md,
            (VRNA_OPTION_MFE | VRNA_OPTION_PF) as u32,
        );
        assert!(!fc.is_null(), "fold_compound returned null");

        let mut mfe_struct: Vec<c_char> = vec![0; n + 1];
        let mfe = vrna_mfe(fc, mfe_struct.as_mut_ptr());

        let mut mfe_scaled: f64 = mfe as f64;
        vrna_exp_params_rescale(fc, &mut mfe_scaled);

        let mut pf_struct: Vec<c_char> = vec![0; n + 1];
        let _pf_energy = vrna_pf(fc, pf_struct.as_mut_ptr());


        let exp_matrices = (*fc).exp_matrices;
        assert!(!exp_matrices.is_null(), "exp_matrices is null");

        let probs_ptr = (*exp_matrices).__bindgen_anon_1.__bindgen_anon_1.probs;
        let iindx_ptr = (*fc).iindx;
        assert!(!probs_ptr.is_null(), "probs is null after vrna_pf");
        assert!(!iindx_ptr.is_null(), "iindx is null");

        let mut defects = vec![0.0_f64; n];
        let mut cost = 0.0_f64;

        for i in 0..n {
            let mut p_paired = 0.0_f64;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let (a, b) = if i < j {
                    (i + 1, j + 1)
                } else {
                    (j + 1, i + 1)
                };
                let iindx_val = *iindx_ptr.add(a) as isize;
                let idx = iindx_val - b as isize;
                let pr_val = *probs_ptr.offset(idx) as f64;
                p_paired += pr_val;
            }
            p_paired = p_paired.min(1.0);

            let should_be_paired = target_bytes[i] != b'.';
            defects[i] = if should_be_paired {
                (1.0 - p_paired).max(0.0)
            } else {
                p_paired
            };
            cost += defects[i];
        }

        let target_c = CString::new(target_no_pk).expect("target has interior NUL");
        let e_target = vrna_eval_structure(fc, target_c.as_ptr()) as f64;

        vrna_fold_compound_free(fc);

        (cost / n as f64, defects, mfe as f64, e_target)
    }
}

fn mutate_ks(seq: &str, pair_map: &HashMap<usize, usize>) -> String {
    let chars: Vec<char> = seq.chars().collect();
    let mut out = chars.clone();
    let mut resolved = vec![false; chars.len()];
    let mut rng = rand::rng();

    for i in 0..chars.len() {
        if resolved[i] {
            continue;
        }
        if !matches!(chars[i], 'K' | 'S') {
            continue;
        }

        let Some(&j) = pair_map.get(&i) else {
            continue;
        };
        if resolved[j] {
            continue;
        }

        match (chars[i], chars[j]) {
            ('K', 'K') => {
                let (a, b) = if rng.random::<f64>() < 0.5 {
                    ('G', 'U')
                } else {
                    ('U', 'G')
                };
                out[i] = a;
                out[j] = b;
            }
            ('S', 'S') => {
                let (a, b) = if rng.random::<f64>() < 0.5 {
                    ('G', 'C')
                } else {
                    ('C', 'G')
                };
                out[i] = a;
                out[j] = b;
            }
            _ => {
                continue;
            }
        }

        resolved[i] = true;
        resolved[j] = true;
    }

    out.into_iter().collect()
}

// pseudoknot-releated functions -> yet to be implemented

fn strip_pseudoknots(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            '[' | ']' => '.',
            _ => c,
        })
        .collect()
}

pub struct PkPrediction {
    pub structure: String,
    pub energy: f64,
    pub dgpk: f64,
    pub dgint: f64,
    pub dg1: f64,
    pub dg2: f64,
    pub start_5: usize,
    pub end_5: usize,
    pub start_3: usize,
    pub end_3: usize,
}

pub fn fold_with_pkplex(seq: &str) -> Vec<PkPrediction> {
    unsafe {
        let mut md: vrna_md_t = std::mem::zeroed();
        vrna_md_set_default(&mut md);
        md.temperature = 37.0;
        md.dangles = 1;

        let seq_c = CString::new(seq).expect("seq has interior NUL");
        let fc = vrna_fold_compound(seq_c.as_ptr(), &md, VRNA_OPTION_MFE as u32);
        assert!(!fc.is_null(), "fold_compound returned null");

        let options = vrna_pk_plex_opt_defaults();
        assert!(
            !options.is_null(),
            "vrna_pk_plex_opt_defaults returned null"
        );

        let accessibility: *mut *const std::os::raw::c_int = std::ptr::null_mut();
        let result_ptr = vrna_pk_plex(fc, accessibility, options);

        let mut out = Vec::new();
        if !result_ptr.is_null() {
            let mut p = result_ptr;
            while !(*p).structure.is_null() {
                let structure = CStr::from_ptr((*p).structure)
                    .to_string_lossy()
                    .into_owned();

                out.push(PkPrediction {
                    structure,
                    energy: (*p).energy,
                    dgpk: (*p).dGpk,
                    dgint: (*p).dGint,
                    dg1: (*p).dG1,
                    dg2: (*p).dG2,
                    start_5: (*p).start_5 as usize,
                    end_5: (*p).end_5 as usize,
                    start_3: (*p).start_3 as usize,
                    end_3: (*p).end_3 as usize,
                });

                libc::free((*p).structure as *mut c_void);
                p = p.add(1);
            }
            libc::free(result_ptr as *mut c_void);
        }

        libc::free(options as *mut c_void);

        vrna_fold_compound_free(fc);
        out
    }
}

// Some additional helper functions

fn find_conserved_positions(sequence: &str, motifs: &[&str]) -> Vec<usize> {
    let mut protected = Vec::new();

    for motif in motifs {
        let motif_len = motif.len();

        if motif_len > sequence.len() {
            continue;
        }

        for start in 0..=sequence.len() - motif_len {
            let window = &sequence[start..start + motif_len];

            if window == *motif {
                for i in start..start + motif_len {
                    protected.push(i);
                }
            }
        }
    }

    protected.sort_unstable();
    protected.dedup();

    protected
}

fn get_local_protected_positions(
    conserved_positions: &[usize],
    slice_start: usize,
    slice_end: usize,
) -> Vec<usize> {
    conserved_positions
        .iter()
        .filter(|&&pos| pos >= slice_start && pos < slice_end)
        .map(|&pos| pos - slice_start)
        .collect()
}

fn assert_slices_balanced(slices: &[Substructure]) -> Result<(), String> {
    for (idx, sub) in slices.iter().enumerate() {
        let opens = sub.structure.matches('(').count();
        let closes = sub.structure.matches(')').count();
        if opens != closes {
            return Err(format!(
                "slice {} range=[{}, {}) is unbalanced: {} '(' vs {} ')'\n  struct: {}",
                idx, sub.start, sub.end, opens, closes, sub.structure
            ));
        }

        let stripped = strip_pseudoknots(&sub.structure);
        let mut depth = 0i32;
        for (pos, c) in stripped.chars().enumerate() {
            match c {
                '(' => depth += 1,
                ')' => {
                    depth -= 1;
                    if depth < 0 {
                        return Err(format!(
                            "slice {} has closing paren before opening at local pos {}",
                            idx, pos
                        ));
                    }
                }
                _ => {}
            }
        }
        if depth != 0 {
            return Err(format!("slice {} not balanced after stripping pk", idx));
        }
    }
    Ok(())
}

fn annotate_pk(seq: &str, target: &str, mfe_struct: &str) -> (String, usize) {
    let pair_map = get_pair_map(target);
    let sb = seq.as_bytes();
    let tb = target.as_bytes();
    let mut s: Vec<char> = mfe_struct.chars().collect();
    let mut n_reinserted = 0;

    for (&i, &j) in pair_map.iter() {
        if tb[i] != b'[' && tb[i] != b']' {
            continue;
        }
        if i >= j {
            continue;
        }

        let ok = matches!(
            (sb[i], sb[j]),
            (b'A', b'U') | (b'U', b'A') | (b'G', b'C') | (b'C', b'G') | (b'G', b'U') | (b'U', b'G')
        );
        if ok {
            if s[i] == '.' {
                s[i] = tb[i] as char;
            }
            if s[j] == '.' {
                s[j] = tb[j] as char;
            }
            if s[i] != '.' && s[j] != '.' {
                n_reinserted += 1;
            }
        }
    }
    (s.into_iter().collect(), n_reinserted)
}

pub fn pk_aware_bp_distance_to_target(seq: &str, target: &str) -> DesignResult {
    let result = bp_distance_to_target(seq, target);
    let pair_map = get_pair_map(target);
    let pk_penalty = pk_pair_mismatches(seq, &pair_map, target) as i64;
    DesignResult {
        bp_distance: result.bp_distance + pk_penalty,
        structure: result.structure,
        mfe: result.mfe,
        sequence: result.sequence,
    }
}

fn identify_mismatches(mfe_structure: &str, target: &str) -> Vec<usize> {
    let target_no_pk = strip_pseudoknots(target);
    let mfe_bytes = mfe_structure.as_bytes();
    let target_bytes = target_no_pk.as_bytes();

    (0..target_bytes.len())
        .filter(|&i| mfe_bytes.get(i).copied() != Some(target_bytes[i]))
        .collect()
}

fn is_adjacent_to_loop(pos: usize, structure: &str) -> bool {
    let bytes = structure.as_bytes();
    let n = bytes.len();
    (pos > 0 && bytes[pos - 1] == b'.') || (pos + 1 < n && bytes[pos + 1] == b'.')
}

fn read_input_file(path: &str) -> Result<(String, String), std::io::Error> {
    let content = fs::read_to_string(path)?;

    let mut seq = String::new();
    let mut target = String::new();

    for line in content.lines() {
        let line = line.trim();

        if line.starts_with("Sequence") {
            if let Some(idx) = line.find(':') {
                seq = line[idx + 1..].trim().to_string();
            }
        } else if line.starts_with("Structure") {
            if let Some(idx) = line.find(':') {
                target = line[idx + 1..].trim().to_string();
            }
        }
    }

    Ok((seq, target))
}

fn insert_ribosome_sequence(seq: &str) -> String {
    //let ribosome_sequence = "GGUUAAGCGACUAAGCGUACACGGUGGAUGCCCUGGCAGUCAGAGGCGAUGAAGGACGUGCUAAUCUGCGAUAAGCGUCGGUAAGGUGAUAUGAACCGUUAUAACCGGCGAUUUCCGAAUGGGGAAACCCAGUGUGUUUCGACACACUAUCAUUAACUGAAUCCAUAGGUUAAUGAGGCGAACCGGGGGAACUGAAACAUCUAAGUACCCCGAGGAAAAGAAAUCAACCGAGAUUCCCCCAGUAGCGGCGAGCGAACGGGGAGCAGCCCAGAGCCUGAAUCAGUGUGUGUGUUAGUGGAAGCGUCUGGAAAGGCGCGCGAUACAGGGUGACAGCCCCGUACACAAAAAUGCACAUGCUGUGAGCUCGAUGAGUAGGGCGGGACACGUGGUAUCCUGUCUGAAUAUGGGGGGACCAUCCUCCAAGGCUAAAUACUCCUGACUGACCGAUAGUGAACCAGUACCGUGAGGGAAAGGCGAAAAGAACCCCGGCGAGGGGAGUGAAAAAGAACCUGAAACCGUGUACGUACAAGCAGUGGGAGCACGCUUAGGCGUGUGACUGCGUACCUUUUGUAUAAUGGGUCAGCGACUUAUAUUCUGUAGCAAGGUUAACCGAAUAGGGGAGCCGAAGGGAAACCGAGUCUUAACUGGGCGUUAAGUUGCAGGGUAUAGACCCGAAACCCGGUGAUCUAGCCAUGGGCAGGUUGAAGGUUGGGUAACACUAACUGGAGGACCGAACCGACUAAUGUUGAAAAAUUAGCGGAUGACUUGUGGCUGGGGGUGAAAGGCCAAUCAAACCGGGAGAUAGCUGGUUCUCCCCGAAAGCUAUUUAGGUAGCGCCUCGUGAAUUCAUCUCCGGGGGUAGAGCACUGUUUCGGCAAGGGGGUCAUCCCGACUUACCAACCCGAUGCAAACUGCGAAUACCGGAGAAUGUUAUCACGGGAGACACACGGCGGGUGCUAACGUCCGUCGUGAAGAGGGAAACAACCCAGACCGCCAGCUAAGGUCCCAAAGUCAUGGUUAAGUGGGAAACGAUGUGGGAAGGCCCAGACAGCCAGGAUGUUGGCUUAGAAGCAGCCAUCAUUUAAAGAAAGCGUAAUAGCUCACUGGUCGAGUCGGCCUGCGCGGAAGAUGUAACGGGGCUAAACCAUGCACCGAAGCUGCGGCAGCGACGCUUAUGCGUUGUUGGGUAGGGGAGCGUUCUGUAAGCCUGCGAAGGUGUGCUGUGAGGCAUGCUGGAGGUAUCAGAAGUGCGAAUGCUGACAUAAGUAACGAUAAAGCGGGUGAAAAGCCCGCUCGCCGGAAGACCAAGGGUUCCUGUCCAACGUUAAUCGGGGCAGGGUGAGUCGACCCCUAAGGCGAGGCCGAAAGGCGUAGUCGAUGGGAAACAGGUUAAUAUUCCUGUACUUGGUGUUACUGCGAAGGGGGGACGGAGAAGGCUAUGUUGGCCGGGCGACGGUUGUCCCGGUUUAAGCGUGUAGGCUGGUUUUCCAGGCAAAUCCGGAAAAUCAAGGCUGAGGCGUGAUGACGAGGCACUACGGUGCUGAAGCAACAAAUGCCCUGCUUCCAGGAAAAGCCUCUAAGCAUCAGGUAACAUCAAAUCGUACCCCAAACCGACACAGGUGGUCAGGUAGAGAAUACCAAGGCGCUUGAGAGAACUCGGGUGAAGGAACUAGGCAAAAUGGUGCCGUAACUUCGGGAGAAGGCACGCUGAUAUGUAGGUGAGGUCCCUCGCGGAUGGAGCUGAAAUCAGUCGAAGAUACCAGCUGGCUGCAACUGUUUAUUAAAAACACAGCACUGUGCAAACACGAAAGUGGACGUAUACGGUGUGACGCCUGCCCGGUGCCGGAAGGUUAAUUGAUGGGGUUAGCGCAAGCGAAGCUCUUGAUCGAAGCCCCGGUAAACGGCGGCCGUAACUAUAACGGUCCUAAGGUAGCGAAAUUCCUUGUCGGGUAAGUUCCGACCUGCACGAAUGGCGUAAUGAUGGCCAGGCUGUCUCCACCCGAGACUCAGUGAAAUUGAACUCGCUGUGAAGAUGCAGUGUACCCGCGGCAAGACGGAAAGACCCCGUGAACCUUUACUAUAGCUUGACACUGAACAUUGAGCCUUGAUGUGUAGGAUAGGUGGGAGGCUUUGAAGUGUGGACGCCAGUCUGCAUGGAGCCGACCUUGAAAUACCACCCUUUAAUGUUUGAUGUUCUAACGUUGACCCGUAAUCCGGGUUGCGGACAGUGUCUGGUGGGUAGUUUGACUGGGGCGGUCUCCUCCUAAAGAGUAACGGAGGAGCACGAAGGUUGGCUAAUCCUGGUCGGACAUCAGGAGGUUAGUGCAAUGGCAUAAGCCAGCUUGACUGCGAGCGUGACGGCGCGAGCAGGUGCGAAAGCAGGUCAUAGUGAUCCGGUGGUUCUGAAUGGAAGGGCCAUCGCUCAACGGAUAAAAGGUACUCCGGGGAUAACAGGCUGAUACCGCCCAAGAGUUCAUAUCGACGGCGGUGUUUGGCACCUCGAUGUCGGCUCAUCACAUCCUGGGGCUGAAGUAGGUCCCAAGGGUAUGGCUGUUCGCCAUUUAAAGUGGUACGCGAGCUGGGUUUAGAACGUCGUGAGACAGUUCGGUCCCUAUCUGCCGUGGGCGCUGGAGAACUGAGGGGGGCUGCUCCUAGUACGAGAGGACCGGAGUGGACGCAUCACUGGUGUUCGGGUUGUCAUGCCAAUGGCACUGCCCGGUAGCUAAAUGCGGAAGAGAUAAGUGCUGAAAGCAUCUAAGCACGAAACUUGCCCCGAGAUGAGUUCUCCCUGACCCUUUAAGGGUCCUGAAGGAACGUUGAAGACGACGACGUUGAUAGGCCGGGUGUGUAAGCGCAGCGAUGCGUUGAGCUAACCGGUACUAAUGAACCGUGAGGCUUAACCU";

    let mut result = seq.as_bytes().to_vec(); // Only use as much of the ribosome sequence as we need 
    let ribosome = &RIBOSOME_SEQUENCE.as_bytes()[..seq.len()];
    for i in 0..result.len() {
        if result[i] == b'N' {
            result[i] = ribosome[i];
        }
    }
    String::from_utf8(result).unwrap()
}

fn mutate_ribosome(seq_in: &str, structure: &str, n_positions: &[usize]) -> String {
    let mut rng = rand::rng();
    let conserved_motifs = important_motifs();
    let pair_map = get_pair_map(structure);
    let conserved_positions = find_conserved_positions(seq_in, &conserved_motifs);
    let hard_loops: bool = find_hard_loops(structure);

    // O(1) membership test for the designable set
    let designable: HashSet<usize> = n_positions.iter().copied().collect();

    let mut comp_dict = HashMap::new();
    comp_dict.insert('A', 'U');
    comp_dict.insert('U', 'A');
    comp_dict.insert('G', 'C');
    comp_dict.insert('C', 'G');

    let _nucleotides = ['A', 'U', 'G', 'C'];
    let paired_nucleotides = ['G', 'C'];
    let purines = ['A', 'G'];

    let mut mut_seq: Vec<char> = seq_in.chars().collect();
    let mut_struct: Vec<char> = structure.chars().collect();

    for i in 0..mut_seq.len() {
        if !designable.contains(&i) {
            continue;
        }

        if conserved_positions.contains(&i) {
            continue;
        }

        if let Some(&j) = pair_map.get(&i) {
            if !designable.contains(&j) {
                // Only fix i from the 5' side to avoid double work
                if i < j {
                    mut_seq[i] = *comp_dict.get(&mut_seq[j]).unwrap();
                }
                continue;
            }

            if i >= j {
                continue;
            }

            if rng.random::<f64>() < 0.20 {
                let choice = paired_nucleotides[rng.random_range(0..paired_nucleotides.len())];
                mut_seq[i] = choice;
                mut_seq[j] = *comp_dict.get(&choice).unwrap();
            }

            if hard_loops && mut_struct[i] == '(' && mut_struct[i + 1] == ')' {
                let choice = if rng.random::<f64>() < 0.5 { 'G' } else { 'C' };
                mut_seq[i] = choice;
                mut_seq[i + 1] = *comp_dict.get(&choice).unwrap();
                println!("HARD LOOP WARNING -> GC-PAIR!");
                if i + 2 < j {
                    if let Some(&partner_of_next) = pair_map.get(&(i + 2)) {
                        if partner_of_next == j - 2
                            && designable.contains(&(i + 2))
                            && designable.contains(&(j - 2))
                        {
                            let choice2 = if rng.random::<f64>() < 0.5 { 'G' } else { 'C' };
                            mut_seq[i + 2] = choice2;
                            mut_seq[j - 2] = *comp_dict.get(&choice2).unwrap();
                        }
                    }
                }
                if i + 3 < j {
                    if let Some(&partner_of_next) = pair_map.get(&(i + 3)) {
                        if partner_of_next == j - 3
                            && designable.contains(&(i + 3))
                            && designable.contains(&(j - 3))
                        {
                            let choice3 = if rng.random::<f64>() < 0.5 { 'G' } else { 'U' };
                            mut_seq[i + 3] = choice3;
                            mut_seq[j - 3] = *comp_dict.get(&choice3).unwrap();
                        }
                    }
                }
            }

            let near_loop = is_adjacent_to_loop(i, structure) || is_adjacent_to_loop(j, structure);

            if near_loop {
                let choice = if rng.random::<f64>() < 0.5 { 'G' } else { 'C' };
                mut_seq[i] = choice;
                mut_seq[j] = *comp_dict.get(&choice).unwrap();

                if i + 1 < j {
                    if let Some(&partner_of_next) = pair_map.get(&(i + 1)) {
                        if partner_of_next == j - 1
                            && designable.contains(&(i + 1))
                            && designable.contains(&(j - 1))
                            && !conserved_positions.contains(&(i + 1))
                            && !conserved_positions.contains(&(j - 1))
                        {
                            let choice2 = if rng.random::<f64>() < 0.5 { 'G' } else { 'C' };
                            mut_seq[i + 1] = choice2;
                            mut_seq[j - 1] = *comp_dict.get(&choice2).unwrap();
                        }
                    }
                }

                if i + 2 < j {
                    if let Some(&partner_of_next) = pair_map.get(&(i + 2)) {
                        if partner_of_next == j - 2
                            && designable.contains(&(i + 2))
                            && designable.contains(&(j - 2))
                            && !conserved_positions.contains(&(i + 2))
                            && !conserved_positions.contains(&(j - 2))
                        {
                            let choice3 = if rng.random::<f64>() < 0.5 { 'G' } else { 'C' };
                            mut_seq[i + 2] = choice3;
                            mut_seq[j - 2] = *comp_dict.get(&choice3).unwrap();
                        }
                    }
                }
            }
        } else {
            // Unpaired designable position
            if rng.random::<f64>() < 0.20 {
                mut_seq[i] = purines[rng.random_range(0..purines.len())];
            }
        }
    }

    mut_seq.into_iter().collect()
}

fn important_motifs() -> &'static [&'static str] {
    let conserved_motifs: &[&str] = &[
        "GGCCAAAU", // Alpha KL
        "UAAACCGG",
        "UUAAUCGAUACCUGGGCUGGCAGAGCGUGCCGGCAUGCUCGGGUGUGAGAUGAGCUGUAUUGAUUGC", // Broccoli
        "CCGUGCGAGACGGUCGGGUCCAUAGCUAAUUCGUUAGUUAUGGAGGCUCGUACGG",             // Broccoli
        "UGAAGCCUCCACG",
        "GCACCUCCGAAGU", // Regular KL
        "AUCACGGGAGACACACGGCGGGUGNNNNNNNNNNNNNNNNNNNNGCAUGCUGGAGGUAUCAGAAGUGCGAAUGCUGACAUAAGUAACGAUAAAGCGGGUGAAAAGCCCGCUCGCCGGAAGACCAAGGGUUCCUGUCCAACGUUAAUCGGGGCAGGGUGAGUCGACCCCUAAGGCGAGGCCGAAAGGCGUAGUCGAUGGGAAACA",
    ];

    return conserved_motifs;
}

fn find_hard_loops(structure: &str) -> bool {
    let motifs: &[&str] = &["()"];
    let mut motif_found: bool = false;
    for n in motifs {
        for i in 0..structure.len() - 1 {
            let j = i + n.len();
            let window = &structure[i..j];
            if window == *n {
                motif_found = true;
                break;
            } else {
                continue;
            };
        }
    }
    motif_found
}

fn ribosome_similarity(seq_in: &str, n_positions: &[usize]) -> RibosomeSimilarity {
    let seq_bytes = seq_in.as_bytes();
    let ref_bytes = RIBOSOME_SEQUENCE.as_bytes();

    assert!(
        seq_bytes.len() <= ref_bytes.len(),
        "seq_in ({} nt) is longer than the reference ribosome sequence ({} nt)",
        seq_bytes.len(),
        ref_bytes.len()
    );

    let mut mismatched_positions = Vec::new();
    for &pos in n_positions {
        assert!(
            pos < seq_bytes.len(),
            "position {} out of bounds for seq_in",
            pos
        );
        if seq_bytes[pos] != ref_bytes[pos] {
            mismatched_positions.push(pos);
        }
    }

    let total = n_positions.len();
    let matches = total - mismatched_positions.len();
    let identity = if total == 0 {
        1.0
    } else {
        matches as f64 / total as f64
    };

    RibosomeSimilarity {
        matches,
        total,
        identity,
        mismatched_positions,
    }
}

fn ask_positive_usize(prompt: &str, default: usize) -> usize {
    loop {
        print!("{prompt} [{default}]: ");
        io::stdout().flush().expect("Failed to flush terminal output");

        let mut input = String::new();
        io::stdin()
            .read_line(&mut input)
            .expect("Failed to read terminal input");

        let input = input.trim();

        // Pressing Enter accepts the default.
        if input.is_empty() {
            return default;
        }

        match input.parse::<usize>() {
            Ok(value) if value > 0 => return value,
            _ => {
                println!("Please enter a positive whole number.");
            }
        }
    }
}

fn ask_positive_i64(prompt: &str, default: i64) -> i64 {
    loop {
        print!("{prompt} [{default}]: ");
        io::stdout().flush().expect("Failed to flush terminal output");

        let mut input = String::new();
        io::stdin()
            .read_line(&mut input)
            .expect("Failed to read terminal input");

        let input = input.trim();

        // Pressing Enter accepts the default.
        if input.is_empty() {
            return default;
        }

        match input.parse::<i64>() {
            Ok(value) if value > 0 => return value,
            _ => {
                println!("Please enter a positive whole number.");
            }
        }
    }
}

fn gc_content(sequence: &str) -> Option<f64> {
    let mut gc_count = 0usize;
    let mut base_count = 0usize;

    for base in sequence.bytes() {
        match base.to_ascii_uppercase() {
            b'G' | b'C' => {
                gc_count += 1;
                base_count += 1;
            }
            b'A' | b'U' => {
                base_count += 1;
            }
            _ => {} 
        }
    }

    if base_count == 0 {
        return None;
    }

    let percentage = (gc_count as f64 / base_count as f64) * 100.0;
    Some((percentage * 100.0).round() / 100.0)
}

fn ask_to_view_results() -> io::Result<bool> {
    print!("\nView all final results in a scrollable window? (y/n): ");
    io::stdout().flush()?;

    let mut answer = String::new();
    io::stdin().read_line(&mut answer)?;

    Ok(answer.trim().eq_ignore_ascii_case("y"))
}

fn show_in_pager(output: &str) -> io::Result<()> {
    let mut pager = Command::new("less")
        .arg("-R")
        .stdin(Stdio::piped())
        .spawn()?;

    if let Some(mut pager_input) = pager.stdin.take() {
        pager_input.write_all(output.as_bytes())?;
    }

    pager.wait()?;
    Ok(())
}

fn gc_cleanup(
    seq_in: &str,
    original_annotation: &str,
    target: &str,
    min_gc_percent: f64,
) -> String {
    assert_eq!(seq_in.len(), original_annotation.len());
    assert_eq!(seq_in.len(), target.len());

    let pair_map = get_pair_map(target);
    let annotation = original_annotation.as_bytes();
    let mut seq: Vec<u8> = seq_in.as_bytes().to_vec();

    let mut current = String::from_utf8(seq.clone()).expect("sequence must be valid UTF-8");
    let mut current_result = bp_distance_to_target(&current, target);
    let mut current_dist = current_result.bp_distance;

    let mut total_replacements = 0usize;
    let mut changed = true;

    while changed {
        changed = false;

        
        let mut pairs: Vec<(usize, usize)> = pair_map
            .iter()
            .filter_map(|(&i, &j)| (i < j).then_some((i, j)))
            .collect();

        pairs.sort_unstable();
        pairs.dedup();

        for (i, j) in pairs {
            
            if annotation[i] != b'N' || annotation[j] != b'N' {
                continue;
            }

            
            let is_gc_pair = matches!(
                (seq[i], seq[j]),
                (b'G', b'C') | (b'C', b'G')
            );

            if !is_gc_pair {
                continue;
            }

            let live_seq = String::from_utf8(seq.clone()).expect("sequence must be valid UTF-8");
            let current_gc = gc_content(&live_seq).unwrap_or(0.0);
            let gc_after = current_gc - (200.0/seq.len() as f64);

            if gc_after < min_gc_percent {
                continue;
            } 

            
            let mut best_trial: Option<(Vec<u8>, DesignResult)> = None;

            for (left, right) in [(b'A', b'U'), (b'U', b'A')] {
                let mut trial = seq.clone();
                trial[i] = left;
                trial[j] = right;

                let trial_string =
                    String::from_utf8(trial.clone()).expect("sequence must be valid UTF-8");
                let trial_result = bp_distance_to_target(&trial_string, target);

                
                if trial_result.bp_distance > current_dist {
                    continue;
                }

                let replace_best = match &best_trial {
                    None => true,
                    Some((_, best_result)) => {
                        trial_result.bp_distance < best_result.bp_distance
                            || (trial_result.bp_distance == best_result.bp_distance
                                && trial_result.mfe < best_result.mfe)
                    }
                };

                if replace_best {
                    best_trial = Some((trial, trial_result));
                }
            }

            if let Some((trial, trial_result)) = best_trial {
                seq = trial;
                current_dist = trial_result.bp_distance;
                current_result = trial_result;
                changed = true;
                total_replacements += 1;
            }
        }
    }

    current = String::from_utf8(seq).expect("sequence must be valid UTF-8");

    println!(
        "GC cleanup complete: {} GC pair(s) converted; final bp_distance={}",
        total_replacements,
        current_result.bp_distance
    );

    current
}

fn ensemble_diversity(seq: &str, target: &str) -> f64 {
    let n = seq.len();
    let target_no_pk = strip_pseudoknots(target);
    let pair_map = get_pair_map(&target_no_pk);

    unsafe {
        let mut md: vrna_md_t = std::mem::zeroed();
        vrna_md_set_default(&mut md);
        md.temperature = 37.0;
        md.dangles = 1;

        let seq_c: CString = CString::new(seq).expect("seq has interior NUL");
        let fc = vrna_fold_compound(seq_c.as_ptr(), &md, VRNA_OPTION_PF as u32);
        assert!(!fc.is_null(), "fold_compound returned null");

        let mut pf_struct: Vec<i8> = vec![0i8; n + 1];
        let _pf_energy = vrna_pf(fc, pf_struct.as_mut_ptr());

        let exp_matrices = (*fc).exp_matrices;
        assert!(!exp_matrices.is_null(), "exp_matrices is null");

        let probs_ptr = (*exp_matrices).__bindgen_anon_1.__bindgen_anon_1.probs;
        let iindx_ptr = (*fc).iindx;
        assert!(!probs_ptr.is_null(), "probs is null after vrna_pf");
        assert!(!iindx_ptr.is_null(), "iindx is null");

        let mut diversity = 0.0_f64;

        for i in 0..n {
            for j in (i + 1)..n {
                let iindx_val: isize = *iindx_ptr.add(i + 1) as isize;
                let idx: isize = iindx_val - (j + 1) as isize;
                let p_ij: f64 = (*probs_ptr.offset(idx)) as f64;

                let is_target_pair: bool = pair_map.get(&i) == Some(&j);

                diversity += if is_target_pair { 1.0 - p_ij } else { p_ij };
            }
        }

        vrna_fold_compound_free(fc);

        diversity
    }
}
