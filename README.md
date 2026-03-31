# tin-nim: Fast Transcript Integrity from BED Coverage

`tin-nim` is a highly efficient command-line tool written in Nim that calculates Transcript Integrity Numbers (TIN) and summary coverage statistics directly from [mosdepth](https://github.com/brentp/mosdepth) per-base coverage files. 

**This tool is NOT a drop-in replacement for RSeQC's [tin.py](https://rseqc.sourceforge.net/#tin-py).** It is designed to be a rapid approximation for pipelines that already generate mosdepth per-base BED files. Because it operates on pre-computed coverage rather than raw BAM alignments, it lacks the read-level context required to perfectly emulate RSeQC. Scores for well-covered transcripts will correlate highly, but sparse or noisy transcripts may diverge significantly.

## Usage

```
Usage:
  main [REQUIRED,optional-params]
Options:
  -h, --help                              print this cligen-erated help
  --help-syntax                           advanced: prepend,plurals,..
  -b=, --bed=           string  REQUIRED  Path to the mosdepth per-base bed.gz file
  -g=, --gtf=           string  REQUIRED  Path to the input GTF annotations file
  -m=, --minCov=        float   1500.0    Minimum total accumulated depth across the transcript to pass
                                          QC (default: 1500.0, approximating 10 reads of 150bp)
  -d=, --dynamicRange=  float   1.0       Log2 ratio between max_cov and min_cov across the transcript
                                          to pass QC (default: 1.0, must be greater than this value to
                                          pass QC)
  --minLength=          int     200       Minimum total length of exons in a given transcript to pass
                                          QC (default: 200, transcripts shorter than this value fail
                                          QC)
  -o=, --output=        string  ""        Optional path to save the TSV output. Defaults to stdout.
```

### Example

```{bash}
tin -b sample.per-base.bed.gz -g gencode.v48.basic.annotation.gtf.gz -o tin.tsv
```

### Installation

Check out the [releases page](https://github.com/jcalendo/tin-nim/releases/) for a pre-compiled linux binary. 

To build the script, a system installation of [htslib](https://www.htslib.org/download/) and of course [Nim](https://nim-lang.org/install.html) are required. Once these dependencies are met:

```{bash}
git clone https://github.com/jcalendo/tin-nim.git
cd tin-nim

# Fetch required Nim packages
nimble install -d 

# Compile with speed optimizations
nim c -d:danger -d:release --opt:speed src/tin.nim
```

## Details

Instead of parsing BAM files and dynamically building pileups, `tin-nim` streams mosdepth's per-base.bed.gz files. It first parses a provided GTF annotation file and builds interval trees using [nim-lapper](https://github.com/brentp/nim-lapper) for every exon in the transcriptome. As it streams the BED file, it queries the trees to find exact overlaps between the continuous coverage blocks and the transcript exons. It then computes TIN using the coverage blocks (i.e. all bases of the exons) using the same core mathematical concept as RSeQC to calculate the effective proportion of the transcript that is uniformly covered: $100 \times e^H / L$. By default, output is streamed as tab-delimited text to stdout but an optional path to a save an output TSV file can also be specified. 

The output TSV file contains one record for each transcript in the supplied GTF. Transcript information includes; `sum_of_exon_lengths` (total number of bases for all exons of the transcript), `fraction_covered` (fraction of bases in the transcript with coverage > 0), `total_cov` (sum of coverage across transcript), `min_cov` (minimum observed coverage across the transcript), `max_cov` (maximum observed coverage across transcript), `mean_cov` (mean coverage across transcript), `median_cov` (median coverage across transcript), `qc_failed` (boolean flag set it transcript did not pass QC filters), `TIN` (transcript integrity number). 

The `qc_failed` flag is currently a minimally tested heuristic that seems to give decently comparable results to RSeQC's tin.py when using only values for which qc_failed == FALSE in the median TIN calculation. **I suggest exploring the resulting output file and playing around with coverage filters especially if you have some tin.py results to compare to**.