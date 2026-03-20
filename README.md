# tin-nim: Fast Transcript Integrity from BED Coverage

`tin-nim` is a highly efficient command-line tool written in Nim that calculates Transcript Integrity Numbers (TIN) and summary coverage statistics directly from [mosdepth](https://github.com/brentp/mosdepth) per-base coverage files. 

**This tool is NOT a drop-in replacement for RSeQC's [tin.py](https://rseqc.sourceforge.net/#tin-py).** It is designed to be a rapid approximation for pipelines that already generate mosdepth BED files. Because it operates on pre-computed coverage  rather than raw BAM alignments, it lacks the read-level context required to perfectly emulate RSeQC. Scores for well-covered transcripts will correlate highly, but sparse or noisy transcripts may diverge significantly.

Instead of parsing BAM files and dynamically building pileups, `tin-nim` streams mosdepth's run-length encoded .bed.gz files.

* It parses a provided GTF annotation file and builds memory-efficient interval trees using [nim-lapper](https://github.com/brentp/nim-lapper) for every exon in the transcriptome. 
* As it streams the BED file, it queries the trees to find exact overlaps between the continuous coverage blocks and the transcript exons.
* It uses the same core mathematical concept as RSeQC to calculate the effective proportion of the transcript that is uniformly covered: $100 \times e^H / L$
* It utilizes the BED file's native run-length encoding to calculate min, max, mean, and median coverage statistics.
