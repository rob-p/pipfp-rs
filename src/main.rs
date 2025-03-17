use bon::Builder;
use clap::Parser;
use niffler::send::from_path;
use pfp::hash::HT;
use pfp::parse::{LT, parse_seq, parse_seq_par};
use rayon::{ThreadPoolBuilder, current_num_threads};
use seq_io::fasta::{self, Record};
use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::time::Instant;
use walkdir::WalkDir;

type Dict = std::collections::HashMap<HT, LT, rustc_hash::FxBuildHasher>;

#[derive(Builder)]
struct CumulativeParse {
    pub dict: Dict,
    pub parses: Vec<HT>,
    pub parse_len: usize,
    pub keep_parses: bool,
}

impl CumulativeParse {
    pub fn new() -> Self {
        Self {
            dict: Dict::default(),
            parses: Vec::<HT>::new(),
            parse_len: 0,
            keep_parses: false,
        }
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (FASTA, possibly compressed)
    #[arg(short, long)]
    input: String,
    /// window size
    #[arg(short, default_value_t = 10)]
    w: usize,
    /// sparsity?
    #[arg(short, default_value_t = 100)]
    p: usize,
    /// Number of threads [default: all]
    #[arg(short, long)]
    threads: Option<usize>,
}

fn merge_parse_in<P: AsRef<Path>>(
    w: usize,
    p: usize,
    path: P,
    global_parse: &mut CumulativeParse,
    threads: usize,
) {
    let (reader, _) = from_path(path).expect("Failed to open input file");
    let mut reader = fasta::Reader::new(reader);
    while let Some(record) = reader.next() {
        let record = record.unwrap();
        let mut seq = Vec::with_capacity(record.seq().len());
        for line in record.seq_lines() {
            seq.extend_from_slice(line);
        }
        let parse = if threads > 1 {
            parse_seq_par(&seq, w, p, threads)
        } else {
            parse_seq(&seq, w, p)
        };

        global_parse.parse_len += parse.phrases.len();
        if global_parse.keep_parses {
            global_parse.parses.extend_from_slice(&parse.phrases);
        }
        global_parse.dict.insert(parse.prefix.0, parse.prefix.1);
        global_parse.dict.insert(parse.suffix.0, parse.suffix.1);
        global_parse
            .dict
            .extend(parse.phrases.iter().copied().zip(parse.phrases_len));
    }
}

fn main() {
    let args = Args::parse();
    let w = args.w;
    let p = args.p;
    let path = Path::new(&args.input);
    let threads = if let Some(t) = args.threads {
        ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .unwrap();
        t
    } else {
        current_num_threads()
    };
    eprintln!(
        "Start prefix-free parsing using {threads} thread{}",
        if threads > 1 { "s" } else { "" }
    );

    let mut cumulative_parse = CumulativeParse::new();
    let mut pi_vals = Vec::<usize>::new();
    let files: Vec<PathBuf> = if path.is_dir() {
        let fa = OsStr::new("fa");
        let fa_cap = OsStr::new("FA");
        WalkDir::new(&path)
            .into_iter()
            .filter_map(Result::ok)
            .filter(|e| match e.path().extension() {
                Some(x) if (x == fa || x == fa_cap) => true,
                _ => false,
            })
            .map(|e| e.path().to_path_buf())
            .collect()
    } else {
        std::fs::read_to_string(&path)
            .unwrap()
            .lines()
            .map(|l| PathBuf::from_str(l).unwrap())
            .collect()
    };

    for entry in files {
        let start_parse = Instant::now();
        merge_parse_in(w, p, entry, &mut cumulative_parse, threads);
        let elapsed = start_parse.elapsed().as_secs_f64();
        eprintln!("Parsed in {:.02} s", elapsed);
        let len_parse = cumulative_parse.parse_len;
        let len_phrases = cumulative_parse
            .dict
            .values()
            .map(|len| *len as usize)
            .sum::<usize>();
        eprintln!("size of dictionary = {}", cumulative_parse.dict.len());
        eprintln!("length of the parse = {len_parse}");
        eprintln!("length of distinct phrases = {len_phrases}");
        eprintln!("pi = {}", len_parse + len_phrases);
        pi_vals.push(len_parse + len_phrases);
    }

    if !pi_vals.is_empty() {
        let tot = *pi_vals.last().expect("non empty") as f64;
        let pi_vals: Vec<f64> = pi_vals.iter().map(|e| (*e as f64) / tot).collect();
        println!("pis = {:?}", pi_vals);
    }
}
