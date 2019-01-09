#include "fasta.hpp"
#include "cluster.hpp"
#include "utils.hpp"
#include "argagg.hpp"
#include "hps/src/hps.h"

#include <iostream>
#include <future>
#include <fstream>
#include <stdio.h>

#include <seqan/align.h>
#include <seqan/graph_msa.h>

// using namespace seqan;

typedef seqan::StringSet<seqan::DnaString> TStringSet;
typedef seqan::Graph<seqan::Alignment<seqan::StringSet<seqan::DnaString, seqan::Dependent<>>, void>> TAlignmentGraph;

int main(int argc, char *argv[]) {
    char* mode = argv[1];
    if (!strcmp(mode, "cluster")) {
        argagg::parser argparser {{
            { "help", {"-h", "--help"},
            "shows this help message", 0},
            { "input", {"-i", "--input"},
            "input fasta/fastq file (required)", 1},
            { "fastq", {"--fastq"},
            "whether input and output should be in fastq format (instead of fasta)", 0},
            { "clusters", {"-c", "--clusters"},
            "clusters file (default: clusters.out)", 1},
            { "threads", {"-t", "--threads"},
            "number of threads to use (default: 1)", 1},
            { "kmer_size", {"-k", "--kmer-size"},
            "k-mer size (default: 6)", 1},
            { "t_s", {"-s", "--score-threshold"},
            "minimum score for two reads to be in the same cluster (default: 0.25)", 1},  
            { "t_v", {"-v", "--max-variance"},
            "max allowed variance for two reads to be in the same cluster (default: 250)", 1},
            { "bv_threshold", {"-B", "--bv-start-threshold"},
            "starting threshold for the bitvector k-mer comparison (default: 0.4)", 1},  
            { "bv_min_threshold", {"-b", "--bv-end-threshold"},
            "ending threshold for the bitvector k-mer comparison (default: 0.2)", 1},  
            { "bv_falloff", {"-f", "--bv-falloff"},
            "falloff value for the bitvector threshold for each iteration (default: 0.05)", 1},  
            { "min_reads_cluster", {"-r", "--min-reads-cluster"},
            "minimum number of reads per cluster (default: 0)", 1},  
        }};

        argagg::parser_results args;
        try {
            args = argparser.parse(argc, argv);
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            return EXIT_FAILURE;
        }

        if (args["help"]) {
            std::cerr << argparser;
            return EXIT_SUCCESS;
        }

        if (!args["input"]) {
            std::cerr << "ERROR: No input file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        std::cerr << "Reading fasta file... ";
        
        // TODO: handle non-existing file
        read_set_t reads;
        if (args["fastq"]) {
            reads = read_fastq_file(args["input"]);
        } else {
            reads = read_fasta_file(args["input"]);
        }

        int n_threads = args["threads"].as<int>(1);
        int kmer_size = args["kmer_size"].as<int>(6);
        double t_s = args["t_s"].as<double>(0.25);
        double t_v = args["t_v"].as<double>(250);
        double bv_threshold = args["bv_threshold"].as<double>(0.4);
        double bv_min_threshold = args["bv_min_threshold"].as<double>(0.2);
        double bv_falloff = args["bv_falloff"].as<double>(0.05);
        int min_reads_cluster = args["min_reads_cluster"].as<int>(0);

        sort_read_set(reads);

        auto clusters = cluster_reads(reads, kmer_size, t_s, t_v, bv_threshold, bv_min_threshold, bv_falloff, min_reads_cluster, n_threads);
        
        std::cerr << "Clustering done" << std::endl;

        // std::cout << clusters.size() << std::endl;
        // int i = 0;
        // for (auto &c : clusters) {
        //     for (auto &cs : c.seqs) {
        //         std::cout << reads[cs.seq_id].header << "," << i << "," << cs.rev << std::endl;
        //     }

        //     ++i;
        // }

        std::ofstream out_file(args["clusters"].as<std::string>("clusters.out"), std::ofstream::binary);
        hps::to_stream(clusters, out_file);
    } else if (!strcmp(mode, "correct")) {
        argagg::parser argparser {{
            { "help", {"-h", "--help"},
            "shows this help message", 0},
            { "input", {"-i", "--input"},
            "input fasta/fastq file (required)", 1},
            { "clusters", {"-c", "--clusters"},
            "clusters file (required)", 1},
            { "fastq", {"--fastq"},
            "whether input and output should be in fastq format (instead of fasta)", 0},
            { "threads", {"-t", "--threads"},
            "number of threads to use (default: 1)", 1},
        }};

        argagg::parser_results args;
        try {
            args = argparser.parse(argc, argv);
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            return EXIT_FAILURE;
        }

        if (args["help"]) {
            std::cerr << argparser;
            return EXIT_SUCCESS;
        }

        if (!args["input"]) {
            std::cerr << "ERROR: No input file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        if (!args["clusters"]) {
            std::cerr << "ERROR: No clusters file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        std::cerr << "Reading fasta file... ";
        
        // TODO: handle non-existing file
        read_set_t reads;
        if (args["fastq"]) {
            reads = read_fastq_file(args["input"]);
        } else {
            reads = read_fasta_file(args["input"]);
        }

        sort_read_set(reads);

        int n_threads = args["threads"].as<int>(1);
        std::ifstream in_file(args["clusters"].as<std::string>(), std::ifstream::binary);
        auto clusters = hps::from_stream<cluster_set_t>(in_file);

        // sort the clusters by size to distribute workload somewhat evenly between cores
        // std::cout << clusters.size() << std::endl;
        // int i = 0;
        // for (auto &c : clusters) {
        //     for (auto &cs : c.seqs) {
        //         std::cout << reads[cs.seq_id].header << "," << i << "," << cs.rev << std::endl;
        //     }

        //     ++i;
        // }

        // msa
        seqan::Score<int> scoringScheme(5, -3, -1, -3); // match, mismatch, extend, open
        auto opts = seqan::MsaOptions<seqan::AminoAcid, seqan::Score<int>>();
        opts.pairwiseAlignmentMethod = 2;
        opts.bandWidth = 500;
        opts.sc = scoringScheme;
        opts.isDefaultPairwiseAlignment = false;
        appendValue(opts.method, 0);
        appendValue(opts.method, 1);

        std::vector<std::future<void>> tasks;
        for (int t = 0; t < n_threads; ++t) {
            tasks.emplace_back(std::async(std::launch::async, [&clusters, n_threads, &reads, opts, t] {
                for (int i = t; i < clusters.size(); i += n_threads) {
                    TStringSet sequences;
                    seqan::StringSet<seqan::String<char>> sequenceNames;

                    for (auto &seq : clusters[i].seqs) {
                        appendValue(sequences, seqan::DnaString(reads[seq.seq_id].seq));
                        appendValue(sequenceNames, seqan::String<char>(reads[seq.seq_id].header));
                    }

                    TAlignmentGraph alignmentGraph(sequences);
                    globalMsaAlignment(alignmentGraph, sequences, sequenceNames, opts);
                    
                    std::cerr << i << "(" << clusters[i].seqs.size() << ") /" << clusters.size() << std::endl;
                }
            }));
        }

        for (auto &&task : tasks) {
            task.get();
        }
        
        // std::cout << alignmentGraph << std::endl;
    } else if (!strcmp(mode, "cluster_summary")) {
        argagg::parser argparser {{
            { "help", {"-h", "--help"},
            "shows this help message", 0},
            { "input", {"-i", "--input"},
            "input fasta/fastq file (required)", 1},
            { "clusters", {"-c", "--clusters"},
            "clusters file (required)", 1},
            { "fastq", {"--fastq"},
            "whether input and output should be in fastq format (instead of fasta)", 0},
        }};

        argagg::parser_results args;
        try {
            args = argparser.parse(argc, argv);
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            return EXIT_FAILURE;
        }

        if (args["help"]) {
            std::cerr << argparser;
            return EXIT_SUCCESS;
        }

        if (!args["input"]) {
            std::cerr << "ERROR: No input file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        if (!args["clusters"]) {
            std::cerr << "ERROR: No clusters file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        std::cerr << "Reading fasta file... ";
        
        // TODO: handle non-existing file
        read_set_t reads;
        if (args["fastq"]) {
            reads = read_fastq_file(args["input"]);
        } else {
            reads = read_fasta_file(args["input"]);
        }

        sort_read_set(reads);

        std::ifstream in_file(args["clusters"].as<std::string>(), std::ifstream::binary);
        auto clusters = hps::from_stream<cluster_set_t>(in_file);

        int cid = 0;
        for (auto c : clusters) {
            for (auto seq : c.seqs) {
                std::cout << reads[seq.seq_id].header << "," << cid << std::endl;
            }

            ++cid;
        }
    } else if (!strcmp(mode, "extract_clusters")) {
        argagg::parser argparser {{
            { "help", {"-h", "--help"},
            "shows this help message", 0},
            { "input", {"-i", "--input"},
            "input fasta/fastq file (required)", 1},
            { "clusters", {"-c", "--clusters"},
            "clusters file (required)", 1},
            { "output", {"-o", "--output-folder"},
            "output folder for fastx files (default: .)", 1},
            { "fastq", {"--fastq"},
            "whether input and output should be in fastq format (instead of fasta)", 0},
        }};

        argagg::parser_results args;
        try {
            args = argparser.parse(argc, argv);
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            return EXIT_FAILURE;
        }

        if (args["help"]) {
            std::cerr << argparser;
            return EXIT_SUCCESS;
        }

        if (!args["input"]) {
            std::cerr << "ERROR: No input file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        if (!args["clusters"]) {
            std::cerr << "ERROR: No clusters file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        std::cerr << "Reading fasta file... ";
        
        // TODO: handle non-existing file
        read_set_t reads;
        if (args["fastq"]) {
            reads = read_fastq_file(args["input"]);
        } else {
            reads = read_fasta_file(args["input"]);
        }

        sort_read_set(reads);

        std::ifstream in_file(args["clusters"].as<std::string>(), std::ifstream::binary);
        auto clusters = hps::from_stream<cluster_set_t>(in_file);

        int cid = 0;
        for (auto c : clusters) {
            std::ostringstream ss_fn;
            if (args["output"]) {
                ss_fn << args["output"].as<std::string>();
                ss_fn << "/";
            }

            ss_fn << "cluster_";
            ss_fn << cid;

            if (args["fastq"]) {
                ss_fn << ".fq";
            } else {
                ss_fn << ".fa";
            }

            std::ofstream cfile;
            cfile.open(ss_fn.str());
            
            for (auto seq : c.seqs) {
                // std::cout << reads[seq.seq_id].header << "," << cid << std::endl;
                cfile << reads[seq.seq_id].header << "\n";
                if (seq.rev) {
                    cfile << reverse_complement(reads[seq.seq_id].seq) << "\n";
                } else {
                    cfile << reads[seq.seq_id].seq << "\n";
                }

                if (args["fastq"]) {
                    cfile << reads[seq.seq_id].ann << "\n";
                    cfile << reads[seq.seq_id].quality << "\n";
                }
            }

            ++cid;
            cfile.close();
        }
    } else {
        std::cerr << "Unknown mode. More info" << std::endl;
    }
    
}
