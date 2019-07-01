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
#include <random>
#include <chrono>
#include <queue>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "Run with mode: ./rattle <cluster|cluster_summary|extract_clusters|correct>" << std::endl;
        return EXIT_FAILURE;
    }

    std::default_random_engine r_eng{static_cast<long unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count())};

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
            "k-mer size for gene clustering (default: 14)", 1},
            { "t_s", {"-s", "--score-threshold"},
            "minimum score for two reads to be in the same gene cluster (default: 0.1)", 1},  
            { "t_v", {"-v", "--max-variance"},
            "max allowed variance for two reads to be in the same gene cluster (default: 500)", 1},
            { "iso", {"--iso"},
            "perform clustering at the isoform level", 0},
            { "iso_kmer_size", {"--iso-kmer-size"},
            "k-mer size for isoform clustering (default: 7)", 1},
            { "iso_t_s", {"--iso-score-threshold"},
            "minimum score for two reads to be in the same isoform cluster (default: 0.25)", 1},  
            { "iso_t_v", {"--iso-max-variance"},
            "max allowed variance for two reads to be in the same isoform cluster (default: 10)", 1},
            { "bv_threshold", {"-B", "--bv-start-threshold"},
            "starting threshold for the bitvector k-mer comparison (default: 0.4)", 1},  
            { "bv_min_threshold", {"-b", "--bv-end-threshold"},
            "ending threshold for the bitvector k-mer comparison (default: 0.2)", 1},  
            { "bv_falloff", {"-f", "--bv-falloff"},
            "falloff value for the bitvector threshold for each iteration (default: 0.05)", 1},  
            { "min_reads_cluster", {"-r", "--min-reads-cluster"},
            "minimum number of reads per cluster (default: 0)", 1},
            { "repr_percentile", {"-p", "--repr-percentile"},
            "cluster representative percentile (default: 0.5)", 1},  
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

        int kmer_size = args["kmer_size"].as<int>(14);
        double t_s = args["t_s"].as<double>(0.1);
        double t_v = args["t_v"].as<double>(500);

        int iso_kmer_size = args["iso_kmer_size"].as<int>(7);
        double iso_t_s = args["iso_t_s"].as<double>(0.25);
        double iso_t_v = args["iso_t_v"].as<double>(10);

        double bv_threshold = args["bv_threshold"].as<double>(0.4);
        double bv_min_threshold = args["bv_min_threshold"].as<double>(0.2);
        double bv_falloff = args["bv_falloff"].as<double>(0.05);

        int min_reads_cluster = args["min_reads_cluster"].as<int>(0);
        double repr_percentile = args["repr_percentile"].as<double>(0.5);

        sort_read_set(reads);
        std::cerr << "Done" << std::endl;

        auto gene_clusters = cluster_reads(reads, kmer_size, t_s, t_v, bv_threshold, bv_min_threshold, bv_falloff, min_reads_cluster, true, repr_percentile, n_threads);
        std::ofstream out_file(args["output"].as<std::string>(".") + "/clusters.out", std::ofstream::binary);
        
        std::cerr << "Gene clustering done" << std::endl;
        std::cerr << gene_clusters.size() << " gene clusters found" << std::endl;
        if (!args["iso"]) {
            hps::to_stream(gene_clusters, out_file);
            return EXIT_SUCCESS;
        }

        // clustering at isoform level
        cluster_set_t iso_clusters;
        int i = 0;
        for (auto &c : gene_clusters) {
            std::cerr << i << " gene cluster to isoform clusters" << std::endl;

            // sort gene cluster seqs by size
            std::stable_sort(c.seqs.begin(), c.seqs.end(), [&reads](cseq_t a, cseq_t b) {
                return reads[a.seq_id].seq.size() > reads[b.seq_id].seq.size();
            });

            // generate new read set with gene cluster reads
            read_set_t gene_reads;
            for (auto &cs : c.seqs) {
                gene_reads.push_back(reads[cs.seq_id]);
            }

            // cluster gene reads & save new iso clusters
            auto iso_clusters_tmp = cluster_reads(gene_reads, iso_kmer_size, iso_t_s, iso_t_v, bv_threshold, bv_min_threshold, bv_falloff, min_reads_cluster, false, repr_percentile, n_threads);
            for (auto &ic : iso_clusters_tmp) {
                cluster_t iso_cluster;
                iso_cluster.main_seq = cseq_t{c.seqs[ic.main_seq.seq_id].seq_id, ic.main_seq.rev};

                for (auto &ics : ic.seqs) {
                    iso_cluster.seqs.push_back(cseq_t{c.seqs[ics.seq_id].seq_id, ics.rev});
                }

                iso_clusters.push_back(iso_cluster);
            }

            ++i;
        }

        hps::to_stream(iso_clusters, out_file);
        return EXIT_SUCCESS;
    } else if (!strcmp(mode, "correct")) {
        argagg::parser argparser {{
            { "help", {"-h", "--help"},
            "shows this help message", 0},
            { "input", {"-i", "--input"},
            "input fasta/fastq file (required)", 1},
            { "clusters", {"-c", "--clusters"},
            "clusters file (required)", 1},
            { "output", {"-o", "--output"},
            "output folder (default: .)", 1},
            { "gap-occ", {"-g", "--gap-occ"},
            "gap-occ (default: 0.6)", 1},
            { "min-occ", {"-m", "--min-occ"},
            "min-occ (default: 0.3)", 1},
            { "split", {"-s", "--split"},
            "split clusters into sub-clusters of size s for msa (default: 100)", 1},
            { "mafft-path", {"--mafft-path"},
            "path to mafft (default: mafft)", 1},
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
        read_set_t reads = read_fastq_file(args["input"]);

        sort_read_set(reads);
        std::cerr << "Done" << std::endl;

        int n_threads = args["threads"].as<int>(1);
        std::ifstream in_file(args["clusters"].as<std::string>(), std::ifstream::binary);
        auto clusters = hps::from_stream<cluster_set_t>(in_file);
        int cid = 0;
        int split = args["split"].as<int>(100);
        double min_occ = args["min-occ"].as<double>(0.3);
        double gap_occ = args["gap-occ"].as<double>(0.6);

        std::queue<pack_to_correct_t> pending_clusters;
        int corrected = 0;
        int total_reads = 0;

        std::ofstream clusters_file;
        clusters_file.open(args["output"].as<std::string>(".") + "/corrected.fq");

        std::ofstream consensi_file;
        consensi_file.open(args["output"].as<std::string>(".") + "/consensi.fa");

        for (auto &tc: clusters) {           
            int n_files = (tc.seqs.size() + split - 1) / split; // ceil(tc.seqs.size / split)

            for (int nf = 0; nf < n_files; nf++) {
                int nreads_in_cluster = (tc.seqs.size() + n_files - 1 - nf) / n_files;
                auto creads = read_set_t(nreads_in_cluster);

                int i = 0;
                for (int j = nf; j < tc.seqs.size(); j += n_files) {
                    auto ts = tc.seqs[j];

                    if (ts.rev) {
                        reads[ts.seq_id].seq = reverse_complement(reads[ts.seq_id].seq);                
                        std::reverse(reads[ts.seq_id].quality.begin(), reads[ts.seq_id].quality.end()); 
                    }

                    creads[i] = reads[ts.seq_id];
                    i++;
                    total_reads++;
                }

                if (creads.size() > 5) {
                    pending_clusters.push(pack_to_correct_t{cid, creads});
                } else {
                    for (int i = 0; i < creads.size(); ++i) {
                        corrected++;
                        clusters_file << creads[i].header << std::endl;
                        clusters_file << creads[i].seq << std::endl;
                        clusters_file << creads[i].ann << std::endl;
                        clusters_file << creads[i].quality << std::endl;
                    }
                }
            }

            cid++;
        }

        std::mutex mu;
        std::vector<std::future<void>> tasks;
        auto consensi = std::vector<read_set_t>(clusters.size());

        for (int t = 0; t < n_threads; ++t) {
            tasks.emplace_back(std::async(std::launch::async, [t, &consensi, &clusters_file, &consensi_file, &pending_clusters, &mu, &corrected, &total_reads, gap_occ, min_occ, n_threads] {
                while (true) {
                    pack_to_correct_t pack;

                    {
                        std::lock_guard<std::mutex> lock(mu);
                        if (pending_clusters.empty()) break;

                        pack = pending_clusters.front();
                        pending_clusters.pop();

                        print_progress(corrected, total_reads);
                    }

                    auto creads = pack.reads;
                    auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(0),
                    5, -4, -8, -6);

                    auto graph = spoa::createGraph();

                    for (int j = 0; j < creads.size(); ++j) {
                        auto alignment = alignment_engine->align(creads[j].seq, graph);
                        graph->add_alignment(alignment, creads[j].seq);
                    }
                    
                    int i = 0;
                    std::vector<std::string> msa;
                    graph->generate_multiple_sequence_alignment(msa);

                    read_set_t aln_reads = read_set_t(creads.size());
                    for (const auto& it: msa) {
                        aln_reads[i] = read_t{creads[i].header, it, "", ""};
                        i++;
                    }

                    auto corrected_reads_pack = correct_reads(creads, aln_reads, min_occ, gap_occ, 30.0, 1);
                    auto corrected_reads = corrected_reads_pack.reads;
                    {
                        std::lock_guard<std::mutex> lock(mu);
                        for (int i = 0; i < corrected_reads.size(); ++i) {
                            clusters_file << corrected_reads[i].header << std::endl;
                            clusters_file << corrected_reads[i].seq << std::endl;
                            clusters_file << corrected_reads[i].ann << std::endl;
                            clusters_file << corrected_reads[i].quality << std::endl;
                        }
                        corrected+=creads.size();
                        consensi[pack.original_cluster_id].push_back(read_t{">r", corrected_reads_pack.consensus, "", ""});
                    }
                }
            }));
        }

        for (auto &&task : tasks) {
            task.get();
        }
        
        print_progress(corrected, total_reads);
        std::cerr << std::endl;

        std::cerr << "Generating consensi..." << std::endl;
        cid = 0;
        for (const auto& it: consensi) {
            if (it.size() > 0) consensi_file << ">cluster_" << cid << " reads=" << it.size() << std::endl;

            if (it.size() > 1) {
                auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(0),
                    5, -4, -8, -6);

                auto graph = spoa::createGraph();

                for (int j = 0; j < it.size(); ++j) {
                    auto alignment = alignment_engine->align(it[j].seq, graph);
                    graph->add_alignment(alignment, it[j].seq);
                }
                
                std::string consensus = graph->generate_consensus();
                consensi_file << consensus << std::endl;
            } else {
                if (it.size() > 0) {
                    consensi_file << it[0].seq << std::endl;
                }
            }

            cid++;
        }

        clusters_file.close();
        consensi_file.close();

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
        
        // TODO: handle non-existing file
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
        
        // TODO: handle non-existing file
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
    } else {
        std::cerr << "Unknown mode. More info" << std::endl;
    }
    
}
