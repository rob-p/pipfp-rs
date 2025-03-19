# `piPFP-rs`: A Rust implementation of piPFP for measuring the growth of sequence collections

`piPFP-rs` is a Rust implementation of `piPFP`[^pipfp]. It is a program to compute a quantity called \(\pi_{w,p}\), which measures the substring complexity of a sequence or 
sequence collection. By looking at how \(\pi_{w,p}\) changes as new sequences are added to the collection, one can infer how much _novel_ sequence is appearing, or how redundant
the newly added sequences are compared to what is already present in the collection.  The parameters \(w\) (the window size) and \(p\) (the sampling factor) affect how  the 
value is computed, but \(\pi\) is fairly robust to these values as described in the preprint of Lipták et al.[^preprint].

## Using `piPFP-rs`

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

## References

[^preprint]: Zsuzsanna Lipták, Simone Lucà, Francesco Masillo bioRxiv 2025.02.21.639270; doi: [https://doi.org/10.1101/2025.02.21.639270](https://doi.org/10.1101/2025.02.21.639270)

[^pipfp]: [pipfp](https://github.com/simolucaa/piPFP)
