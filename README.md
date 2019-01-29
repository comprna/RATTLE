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

## Clustering
```
$ ./rattle cluster
ERROR: No input file provided
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
        k-mer size (default: 6)
    -s, --score-threshold
        minimum score for two reads to be in the same cluster (default: 0.25)
    -v, --max-variance
        max allowed variance for two reads to be in the same cluster (default: 250)
    -B, --bv-start-threshold
        starting threshold for the bitvector k-mer comparison (default: 0.4)
    -b, --bv-end-threshold
        ending threshold for the bitvector k-mer comparison (default: 0.2)
    -f, --bv-falloff
        falloff value for the bitvector threshold for each iteration (default: 0.05)
    -r, --min-reads-cluster
        minimum number of reads per cluster (default: 0)
```

Recommended parameters for ONT data clustering at the gene level are `-k 14 -s 0.1 -v 500` and for clustering at the transcript level `-k 7 -s 0.25 -v 10`

This clustering step will generate a file containing its clusters in binary format (by default clusters.out). To work with these clusters, the following two commands are used:

## Clustering summary
```
$ ./rattle cluster_summary
ERROR: No input file provided
    -h, --help
        shows this help message
    -i, --input
        input fasta/fastq file (required)
    -c, --clusters
        clusters file (required)
    --fastq
        whether input and output should be in fastq format (instead of fasta)
```

This command will generate a CSV file containing two columns:
```
read_id, cluster_id
```

## Clustering extraction
```
$ ./rattle extract_clusters
ERROR: No input file provided
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
$ ./rattle correct
ERROR: No input file provided
    -h, --help
        shows this help message
    -i, --input
        input fasta/fastq file (required)
    -c, --clusters
        clusters file (required)
    --fastq
        whether input and output should be in fastq format (instead of fasta)
    -t, --threads
        number of threads to use (default: 1)
```

This command will perform error correction based on a clustering. This is a WIP and is still not working as of 09/01/2019
