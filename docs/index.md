# `pipfp-rs`: A Rust implementation of piPFP for measuring the growth of sequence collections

`pipfp-rs` is a Rust implementation of `pipfp`[^pipfp]. It is a program to compute a quantity called \(\pi_{w,p}\), which measures the substring complexity of a sequence or 
sequence collection. By looking at how \(\pi_{w,p}\) changes as new sequences are added to the collection, one can infer how much _novel_ sequence is appearing, or how redundant
the newly added sequences are compared to what is already present in the collection.  The parameters \(w\) (the window size) and \(p\) (the sampling factor) affect how  the 
value is computed, but \(\pi\) is fairly robust to these values as described in the preprint of Lipták et al.[^preprint].

## How `pipfp-rs` works

The usage for `pipfp-rs` is shown below.  It has one required parameter (`-i`/`--input`) which takes either the path to a directory, or a file that lists a collection of input files (1 file per line).
If it is provided with a directory, it will walk the directory and collect all of the files ending in `.fa`, `.fa.gz`, `.fna` or `.fna.gz` (or the capital versions of these), it will then sort the 
file names lexicographically, and process the sequences in that order. If a list of files is provided, the files will be processed in the order they are given.

To process the input `pipfp-rs` will iterate through the files, and after parsing each, will compute the current value of \(\pi_{w,p}\) --- the substring complexity measure for the current
sequence collection.  These \(\pi\) values will be stored and, at the end of the program, written as a `JSON` dictinary to `stdout`.

## Using `pipfp-rs`

The usage of `pipfp-rs` is as follows.

```sh
Usage: pipfp-rs [OPTIONS] --input <INPUT>

Options:
  -i, --input <INPUT>      Either a directory containing FASTA files, or a file with a list of (FASTA) file paths
  -w <W>                   window size [default: 10]
  -p <P>                   sparsity [default: 100]
  -t, --threads <THREADS>  Number of threads [default: all]
  -n, --normalized         write out normalized π values
  -h, --help               Print help
  -V, --version            Print version
```

The `-i` parameter is the only required input.  The `-w` and `-p` parameters control the window size and sparsity factor, as describe in the pipfp preprint[^preprint]. The 
`-t` parameter allows for multi-threaded comptuation of the \(pi\) values, and can speed up execution if the input sequences are large or there are many.  Finally, the `-n`
flag will write out the \(pi\) values in normalized form (i.e. so that the largest value at the end of the array is 1.0).

The output will be written to `stdout`, and consists of a `JSON` dictionary with a key `pis` whose corresponding value is an array of the compute \(\pi\) values and 
a key `norm` whose value is `true` if the `-n` flag was used and `false` otherwise.

## References

[^preprint]: Zsuzsanna Lipták, Simone Lucà, Francesco Masillo bioRxiv 2025.02.21.639270; doi: [https://doi.org/10.1101/2025.02.21.639270](https://doi.org/10.1101/2025.02.21.639270)

[^pipfp]: [pipfp](https://github.com/simolucaa/piPFP)
