#include "fasta.hpp"
#include "cluster.hpp"
#include "utils.hpp"
#include "correct.hpp"
#include "argagg.hpp"
#include "hps/src/hps.h"
#include "spoa/spoa.hpp"

#include <iostream>
#include <future>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <queue>
#include <unistd.h>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "Run with mode: ./rattle <cluster|cluster_summary|extract_clusters|correct|polish>" << std::endl;
        return EXIT_FAILURE;
    }

    char* mode = argv[1];
    if (!strcmp(mode, "cluster")) {
        argagg::parser argparser {{
            { "help", {"-h", "--help"},
            "shows this help message", 0},
            { "input", {"-i", "--input"},
            "input fasta/fastq file (required)", 1},
            { "fastq", {"--fastq"},
            "whether input and output should be in fastq format (instead of fasta)", 0},
            { "output", {"-o", "--output"},
            "output folder (default: .)", 1},
            { "threads", {"-t", "--threads"},
            "number of threads to use (default: 1)", 1},
            { "kmer_size", {"-k", "--kmer-size"},
            "k-mer size for gene clustering (default: 10)", 1},
            { "t_s", {"-s", "--score-threshold"},
            "minimum score for two reads to be in the same gene cluster (default: 0.2)", 1},  
            { "t_v", {"-v", "--max-variance"},
            "max allowed variance for two reads to be in the same gene cluster (default: 1000000)", 1},
            { "iso", {"--iso"},
            "perform clustering at the isoform level", 0},
            { "iso_kmer_size", {"--iso-kmer-size"},
            "k-mer size for isoform clustering (default: 11)", 1},
            { "iso_t_s", {"--iso-score-threshold"},
            "minimum score for two reads to be in the same isoform cluster (default: 0.3)", 1},  
            { "iso_t_v", {"--iso-max-variance"},
            "max allowed variance for two reads to be in the same isoform cluster (default: 25)", 1},
            { "bv_threshold", {"-B", "--bv-start-threshold"},
            "starting threshold for the bitvector k-mer comparison (default: 0.4)", 1},  
            { "bv_min_threshold", {"-b", "--bv-end-threshold"},
            "ending threshold for the bitvector k-mer comparison (default: 0.2)", 1},  
            { "bv_falloff", {"-f", "--bv-falloff"},
            "falloff value for the bitvector threshold for each iteration (default: 0.05)", 1},  
            { "min_reads_cluster", {"-r", "--min-reads-cluster"},
            "minimum number of reads per cluster (default: 0)", 1},
            { "repr_percentile", {"-p", "--repr-percentile"},
            "cluster representative percentile (default: 0.15)", 1},  
            { "rna", {"--rna"},
            "use this mode if data is direct RNA (disables checking both strands)", 0},
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

        int n_threads = args["threads"].as<int>(1);

        int kmer_size = args["kmer_size"].as<int>(10);
        double t_s = args["t_s"].as<double>(0.2);
        double t_v = args["t_v"].as<double>(1000000);

        int iso_kmer_size = args["iso_kmer_size"].as<int>(11);
        double iso_t_s = args["iso_t_s"].as<double>(0.3);
        double iso_t_v = args["iso_t_v"].as<double>(25);

        double bv_threshold = args["bv_threshold"].as<double>(0.4);
        double bv_min_threshold = args["bv_min_threshold"].as<double>(0.2);
        double bv_falloff = args["bv_falloff"].as<double>(0.05);

        int min_reads_cluster = args["min_reads_cluster"].as<int>(0);
        double repr_percentile = args["repr_percentile"].as<double>(0.15);

        if(kmer_size > 16 || iso_kmer_size > 16){
            std::cerr << "\nError: maximum kmer size = 16 \n";
            return EXIT_FAILURE;
        }

        bool is_rna = args["rna"];
        std::cerr << "RNA mode: " << std::boolalpha << is_rna << std::endl;

        std::cerr << "Reading fasta file... ";
        read_set_t reads;
        if(access(args["input"], F_OK )){
            std::cerr << "\nError: Input file not found! \n";
            return EXIT_FAILURE;
        } else {
            std::string filename = args["input"];
            int i = filename.find_last_of(".");
            std::string extension = filename.substr(i + 1);
            if (!extension.compare("fq") || !extension.compare("fastq")){
                reads = read_fastq_file(args["input"]);
            } else if (!extension.compare("fasta") || !extension.compare("fa")){
                reads = read_fasta_file(args["input"]);
            } else {
                std::cerr << "\nError: Input file format incorrect! Please use fasta/fastq file. \n";
                return EXIT_FAILURE;
            }
        }

        sort_read_set(reads);
        std::cerr << "Done" << std::endl;

        auto gene_clusters = cluster_reads(reads, kmer_size, t_s, t_v, bv_threshold, bv_min_threshold, bv_falloff, min_reads_cluster, false, repr_percentile, is_rna, true, n_threads);
        std::ofstream out_file(args["output"].as<std::string>(".") + "/clusters.out", std::ofstream::binary);
        
        std::cerr << "Gene clustering done" << std::endl;
        std::cerr << gene_clusters.size() << " gene clusters found" << std::endl;
        if (!args["iso"]) {
            hps::to_stream(gene_clusters, out_file);
            out_file.close();
            return EXIT_SUCCESS;
        }

        // clustering at isoform level
        cluster_set_t iso_clusters;
        int i = 0;
        for (auto &c : gene_clusters) {
            // sort gene cluster seqs by size
            std::stable_sort(c.seqs.begin(), c.seqs.end(), [&reads](cseq_t a, cseq_t b) {
                return a.seq_id > b.seq_id;
            });

            std::stable_sort(c.seqs.begin(), c.seqs.end(), [&reads](cseq_t a, cseq_t b) {
                return reads[a.seq_id].seq.size() > reads[b.seq_id].seq.size();
            });

            // generate new read set with gene cluster reads
            read_set_t gene_reads;
            for (auto &cs : c.seqs) {
                gene_reads.push_back(reads[cs.seq_id]);
            }

            // cluster gene reads & save new iso clusters
            auto iso_clusters_tmp = cluster_reads(gene_reads, iso_kmer_size, iso_t_s, iso_t_v, bv_threshold, bv_min_threshold, bv_falloff, min_reads_cluster, false, repr_percentile, is_rna, false, n_threads);
            for (auto &ic : iso_clusters_tmp) {
                cluster_t iso_cluster;
                iso_cluster.main_seq = cseq_t{c.seqs[ic.main_seq.seq_id].seq_id, ic.main_seq.rev};

                for (auto &ics : ic.seqs) {
                    iso_cluster.seqs.push_back(cseq_t{c.seqs[ics.seq_id].seq_id, ics.rev});
                }

                iso_clusters.push_back(iso_cluster);
            }

            ++i;
            print_progress(i, gene_clusters.size());
        }

        std::cerr << "Isoform clustering done" << std::endl;
        std::cerr << iso_clusters.size() << " isoform clusters found" << std::endl;
        hps::to_stream(iso_clusters, out_file);
        out_file.close();
        return EXIT_SUCCESS;
    } else if (!strcmp(mode, "correct")) {
        argagg::parser argparser {{
            { "help", {"-h", "--help"},
            "shows this help message", 0},
            { "input", {"-i", "--input"},
            "input fasta/fastq file (required)", 1},
            { "fastq", {"--fastq"},
            "whether input and output should be in fastq format (instead of fasta)", 0},
            { "clusters", {"-c", "--clusters"},
            "clusters file (required)", 1},
            { "output", {"-o", "--output"},
            "output folder (default: .)", 1},
            { "gap-occ", {"-g", "--gap-occ"},
            "gap-occ (default: 0.3)", 1},
            { "min-occ", {"-m", "--min-occ"},
            "min-occ (default: 0.3)", 1},
            { "split", {"-s", "--split"},
            "split clusters into sub-clusters of size s for msa (default: 200)", 1},
            { "min-reads", {"-r", "--min-reads"},
            "min reads to correct/output consensus for a cluster (default: 5)", 1},
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
        read_set_t reads;
        if(access(args["input"], F_OK )){
            std::cerr << "\nError: Input file not found! \n";
            return EXIT_FAILURE;
        } else {
            std::string filename = args["input"];
            int i = filename.find_last_of(".");
            std::string extension = filename.substr(i + 1);
            if (!extension.compare("fq") || !extension.compare("fastq")){
                reads = read_fastq_file(args["input"]);
            } else if (!extension.compare("fasta") || extension.compare("fa")){
                reads = read_fasta_file(args["input"]);
            } else {
                std::cerr << "\nError: Input file format incorrect! Please use fasta/fastq file. \n";
                return EXIT_FAILURE;
            }
        }

        sort_read_set(reads);
        std::cerr << "Done" << std::endl;

        int n_threads = args["threads"].as<int>(1);
        std::ifstream in_file(args["clusters"].as<std::string>(), std::ifstream::binary);
        auto clusters = hps::from_stream<cluster_set_t>(in_file);
        int split = args["split"].as<int>(200);
        double min_occ = args["min-occ"].as<double>(0.3);
        double gap_occ = args["gap-occ"].as<double>(0.3);
        int min_reads = args["min-reads"].as<int>(5);

        correction_results_t correction = correct_reads(clusters, reads, min_occ, gap_occ, 30.0, split, min_reads, n_threads);
        write_fastq_file(correction.corrected, args["output"].as<std::string>(".") + "/corrected.fq");
        write_fastq_file(correction.uncorrected, args["output"].as<std::string>(".") + "/uncorrected.fq");
        write_fastq_file(correction.consensi, args["output"].as<std::string>(".") + "/consensi.fq");

        std::cerr << "Done" << std::endl;
        
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
        if(access(args["input"], F_OK )){
            std::cerr << "\nError: Input file not found! \n";
            return EXIT_FAILURE;
        }
        
        read_set_t reads;
        if (args["fastq"]) {
            reads = read_fastq_file(args["input"]);
        } else {
            reads = read_fasta_file(args["input"]);
        }

        sort_read_set(reads);
        std::cerr << "Done" << std::endl;

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
            { "minreads", {"-m", "--min-reads"},
            "min reads per cluster to save it into a file", 1},
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
        if(access(args["input"], F_OK )){
            std::cerr << "\nError: Input file not found! \n";
            return EXIT_FAILURE;
        }
        
        read_set_t reads;
        if (args["fastq"]) {
            reads = read_fastq_file(args["input"]);
        } else {
            reads = read_fasta_file(args["input"]);
        }

        sort_read_set(reads);
        std::cerr << "Done" << std::endl;

        std::ifstream in_file(args["clusters"].as<std::string>(), std::ifstream::binary);
        auto clusters = hps::from_stream<cluster_set_t>(in_file);
        int min_reads = args["minreads"].as<int>(0);

        int cid = 0;
        for (auto c : clusters) {
            if (c.seqs.size() > min_reads) {
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
                
                cfile.close();
            }

            ++cid;
        }
    } else if (!strcmp(mode, "polish")) {
        argagg::parser argparser {{
            { "help", {"-h", "--help"},
            "shows this help message", 0},
            { "input", {"-i", "--input"},
            "input RATTLE consensi fasta/fastq file (required)", 1},
            { "output", {"-o", "--output-folder"},
            "output folder for fastx files (default: .)", 1},
            { "threads", {"-t", "--threads"},
            "number of threads to use (default: 1)", 1},
            { "rna", {"--rna"},
            "use this mode if data is direct RNA (disables checking both strands)", 0},
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
        if(access(args["input"], F_OK )){
            std::cerr << "\nError: Input file not found! \n";
            return EXIT_FAILURE;
        }
        
        read_set_t reads = read_fastq_file(args["input"]);

        sort_read_set(reads);
        std::cerr << "Done" << std::endl;

        int n_threads = args["threads"].as<int>(1);
        bool is_rna = args["rna"];

        std::cerr << "Clustering consensus sequences..." << std::endl;
        auto clusters = cluster_reads(reads, 6, 0.5, 25, 0.4, 0.4, 0.05, 0, false, 0.15, is_rna, true, n_threads);
        auto correction = correct_reads(clusters, reads, 0.3, 0.3, 30.0, 200, 0, n_threads);

        int cid = 0;
        for (auto &r: correction.consensi) {
            int total_reads = 0;
            auto creads = clusters[cid].seqs;

            for (auto &s: creads) {
                auto info = split(reads[s.seq_id].header, '=');
                int rcount = std::stoi(info[1]);
                total_reads += rcount;
            }

            r.header += " total_reads=" + std::to_string(total_reads);
            cid++;
        }

        write_fastq_file(correction.consensi, args["output"].as<std::string>(".") + "/transcriptome.fq");
    } else {
        std::cerr << "Unknown mode. More info" << std::endl;
    }
    
}