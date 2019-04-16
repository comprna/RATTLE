# RATTLE
Reference-free reconstruction and error correction of transcriptomes from Nanopore long-read sequencing

# Requirements
GCC with C++14 suppport

# Installation
* Clone the repository
```
git clone --recurse-submodules https://github.com/comprna/RATTLE
```

* Build RATTLE
```
cd RATTLE
make
```

# Running RATTLE
**Warning** All the commands and parameters are still highly experimental and subject to changes in future versions

## Sample pipeline
Details on each command can be found below, but these are the most common commands used when running RATTLE.

* Cluster cDNA/RNA Nanopore reads at gene level
```
$ ./rattle cluster -i reads.fq -c genes.out --fastq
```

* Cluster cDNA/RNA Nanopore reads at isoform level
```
$ ./rattle cluster -i reads.fq -c transcripts.out --fastq --iso
```

* View clustering summary (csv with read_id,cluster_id)
```
$ ./rattle cluster_summary -i reads.fq -c transcripts.out --fastq
```

* Extract 1 fastq file per cluster in clusters/ folder
```
$ mkdir clusters
$ ./rattle extract_clusters -i reads.fq -c transcripts.out -o clusters --fastq
```

* Correct reads (mafft in $PATH)
```
$ ./rattle correct -i reads.fq -c transcripts.out > corrected.fq
```

* Extract 1 fastq file per corrected cluster in corrected_clusters/ folder
```
$ mkdir corrected_clusters
$ ./rattle extract_clusters -i corrected.fq -c transcripts.out -o corrected_clusters --fastq
```


## Clustering
```
$ ./rattle cluster -h
    -h, --help
        shows this help message
    -i, --input
        input fasta/fastq file (required)
    --fastq
        whether input and output should be in fastq format (instead of fasta)
    -c, --clusters
        clusters file (default: clusters.out)
    -t, --threads
        number of threads to use (default: 1)
    -k, --kmer-size
        k-mer size for gene clustering (default: 14)
    -s, --score-threshold
        minimum score for two reads to be in the same gene cluster (default: 0.1)
    -v, --max-variance
        max allowed variance for two reads to be in the same gene cluster (default: 500)
    --iso
        perform clustering at the isoform level
    --iso-kmer-size
        k-mer size for isoform clustering (default: 7)
    --iso-score-threshold
        minimum score for two reads to be in the same isoform cluster (default: 0.25)
    --iso-max-variance
        max allowed variance for two reads to be in the same isoform cluster (default: 10)
    -B, --bv-start-threshold
        starting threshold for the bitvector k-mer comparison (default: 0.4)
    -b, --bv-end-threshold
        ending threshold for the bitvector k-mer comparison (default: 0.2)
    -f, --bv-falloff
        falloff value for the bitvector threshold for each iteration (default: 0.05)
    -r, --min-reads-cluster
        minimum number of reads per cluster (default: 0)
```

This clustering step will generate a file containing its clusters in binary format (by default clusters.out). To work with these clusters, the following commands are used.

## Clustering summary
```
$ ./rattle cluster_summary -h
    -h, --help
        shows this help message
    -i, --input
        input fasta/fastq file (required)
    -c, --clusters
        clusters file (required)
    --fastq
        whether input and output should be in fastq format (instead of fasta)

```

This command will output a CSV file to stdout containing two columns:
```
read_id, cluster_id
```

## Clustering extraction
```
$ ./rattle extract_clusters -h
    -h, --help
        shows this help message
    -i, --input
        input fasta/fastq file (required)
    -c, --clusters
        clusters file (required)
    -o, --output-folder
        output folder for fastx files (default: .)
    --fastq
        whether input and output should be in fastq format (instead of fasta)

```

This command will generate 1 fastq/fasta file per cluster in the specified folder. This is useful to analyze particular clusters, perform alignments...

## Error correction
```
$ ./rattle correct -h
    -h, --help
        shows this help message
    -i, --input
        input fastq file (required)
    -c, --clusters
        clusters file (required)
    --mafft-path
        path to mafft (default: mafft)
    -t, --threads
        number of threads to use (default: 1)
```

This command will perform error correction based on a previous clustering using MSA from MAFFT [1]. Path to mafft can be specified using the `--mafft-path` argument if it is not in `$PATH`. Input must be fastq (quality values are used in the correction step).
