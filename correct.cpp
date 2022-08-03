#include <vector>
#include <future>
#include <mutex>
#include <string>
#include <algorithm>
#include <iostream>
#include <queue>
#include <fstream>

#include "correct.hpp"
#include "cluster.hpp"
#include "utils.hpp"
#include "spoa/spoa.hpp"

void print_vector(const std::vector<char> &v) {
    for (int i = 0; i < v.size(); ++i) std::cout<<int(v[i]) << " ";
    std::cout << std::endl;
}

std::vector<std::string> splitString(std::string str, char delimiter) {
    std::vector<std::string> internal;
    std::stringstream ss(str); // Turn the string into a stream.
    std::string tok;

    while(getline(ss, tok, delimiter)) {
        internal.push_back(tok);
    }

    return internal;
}

void fix_msa_ends(read_set_t &reads, msa_t &aln) {
    for (int i = 0; i < aln.size(); ++i) {
        bool reversed = false;

remove_blocks:
        int pos = 0;
        int end_pos = 0;
        while (pos < aln[i].size()) {
            while (pos < aln[i].size() && aln[i][pos] == '-') ++pos;
            
            end_pos = pos;
            int gaps = 0;
            int sz = 0;
            while (gaps < 4 && end_pos < aln[i].size()) {
                if (aln[i][end_pos] == '-') ++gaps;
                else {
                    ++sz;
                    gaps = 0;
                }

                ++end_pos;
            }

            if (sz < 10) {
                // find large gap after small block
                while (end_pos < aln[i].size() && aln[i][end_pos] == '-') {
                    ++end_pos;
                    ++gaps;
                } 
                
                if (gaps >= 20) {
                    for (int j = pos; j < end_pos; ++j) {
                        aln[i][j] = '-';
                    }

                    reads[i].quality.erase(0, sz);
                    reads[i].seq.erase(0, sz);
                    pos = end_pos;
                } else {
                    std::reverse(aln[i].begin(), aln[i].end());
                    std::reverse(reads[i].quality.begin(), reads[i].quality.end());
                    std::reverse(reads[i].seq.begin(), reads[i].seq.end());
                    if (!reversed) {
                        reversed = true;
                        goto remove_blocks;
                    }
                    break;
                }
            } else {
                std::reverse(aln[i].begin(), aln[i].end());
                std::reverse(reads[i].quality.begin(), reads[i].quality.end());
                std::reverse(reads[i].seq.begin(), reads[i].seq.end());
                if (!reversed) {
                    reversed = true;
                    goto remove_blocks;
                }
                break;
            }
        }
    }
}

consensus_vector_t generate_consensus_vector(const read_set_t &reads, const msa_t &aln, int n_threads) {
    // Empty reads check
    if(reads.size() == 0 || aln.size() == 0){
        return {};
    }

    // generate consensus vector
    auto nt_info = std::vector<map_nt_info_t>(aln[0].size());
    for (int i = 0; i < nt_info.size(); i++) {
    // for (int i = nt_info.size() - 1; i >= 0; i--) {
        nt_info[i] = map_nt_info_t();
        nt_info[i]['A'] = pos_info_t{'A', 0.0, 0, 0};
        nt_info[i]['C'] = pos_info_t{'C', 0.0, 0, 0};
        nt_info[i]['T'] = pos_info_t{'T', 0.0, 0, 0};
        nt_info[i]['U'] = pos_info_t{'U', 0.0, 0, 0};
        nt_info[i]['G'] = pos_info_t{'G', 0.0, 0, 0};
        nt_info[i]['-'] = pos_info_t{'-', 0.0, 0, 0};
    }

    std::mutex mu;

    std::vector<std::future<void>> tasks;
    for (int t = 0; t < n_threads; ++t) {
        tasks.emplace_back(std::async(std::launch::async, [t, &reads, &aln, n_threads, &mu, &nt_info] {
            std::vector<map_nt_info_t> local_nt_info = std::vector<map_nt_info_t>(aln[0].size());
            for (int i = 0; i < local_nt_info.size(); i++) {
                local_nt_info[i] = map_nt_info_t();
                local_nt_info[i]['A'] = pos_info_t{'A', 0.0, 0, 0};
                local_nt_info[i]['C'] = pos_info_t{'C', 0.0, 0, 0};
                local_nt_info[i]['T'] = pos_info_t{'T', 0.0, 0, 0};
                local_nt_info[i]['U'] = pos_info_t{'U', 0.0, 0, 0};
                local_nt_info[i]['G'] = pos_info_t{'G', 0.0, 0, 0};
                local_nt_info[i]['-'] = pos_info_t{'-', 0.0, 0, 0};
            }
            // std::vector<map_nt_info_t> local_nt_info = nt_info;

            for (int i = t; i < reads.size(); i+=n_threads) {
                auto read_aln = aln[i];
                int seq_pos = -1;

                for (int k = 0; k < read_aln.size(); k++) {
                    char nt = read_aln[k];
                    double err_p = 0.0;
                    if (nt != '-') {
                        seq_pos++;
                        err_p = phred_err(reads[i].quality[seq_pos]);
                    }

                    if (seq_pos >= 0 && seq_pos < reads[i].quality.size()) {
                        local_nt_info[k][nt].occ++;
                        local_nt_info[k][nt].err += err_p;

                        if (seq_pos == reads[i].quality.size() - 1) {
                            seq_pos++; // end of read
                            // break;
                        }
                    }
                }
            }

            std::lock_guard<std::mutex> lock(mu);
            for (int k = 0; k < aln[0].size(); ++k) {
                for (auto& kv : local_nt_info[k]) {
                    nt_info[k][kv.first].occ += kv.second.occ;
                    nt_info[k][kv.first].err += kv.second.err;
                }
            }
        }));
    }

    for (auto &&task : tasks) {
        task.get();
    }

    // generate mean error and consensus vector
    // std::cerr << "Aln0 ss: " << aln[0].seq.size() << std::endl;
    auto consensus_nt = std::vector<char>(aln[0].size());
    for (int k = 0; k < aln[0].size(); ++k) {
        int max_occ = 0;
        char max_nt = 0;
        for (auto& kv : nt_info[k]) {
            if (kv.second.occ > 0) {
                for (auto& kv2 : nt_info[k]) {
                    nt_info[k][kv.first].total_occ += kv2.second.occ;
                }

                nt_info[k][kv.first].err /= double(kv.second.occ);
            }

            if (kv.second.occ > max_occ) {
                max_occ = kv.second.occ;
                max_nt = kv.first;
            }
        }

        if (max_nt == 0) max_nt = '-';
        consensus_nt[k] = max_nt;
    }
    return consensus_vector_t{nt_info, consensus_nt};
}

// TODO: n_threads is useless here, default set to 1
corrected_pack_t correct_read_pack(const read_set_t &reads, const msa_t &aln, double min_occ, double gap_occ, double err_ratio, int n_threads) {
    auto cv = generate_consensus_vector(reads, aln, n_threads);
    auto nt_info = cv.nt_info;
    auto consensus_nt = cv.consensus_nt;

    // correct aln
    read_set_t corrected_reads;
    read_set_t uncorrected_reads;
    std::mutex mu;
    std::vector<std::future<void>> tasks;

    for (int t = 0; t < n_threads; ++t) {
        // tasks.emplace_back(std::async(std::launch::async, [t, &reads, &aln, n_threads, &mu, &nt_info, &consensus_nt, min_occ, gap_occ, err_ratio, &corrected_reads, &uncorrected_reads] {
            for (int i = t; i < reads.size(); i+=n_threads) {
                auto read_aln = aln[i];
                int seq_pos = -1;
                std::string res_read;
                std::string res_qt;

                int n2g = 0;
                int g2n = 0;
                int n2n = 0;

                for (int k = 0; k < read_aln.size(); ++k) {
                    char nt = read_aln[k];
                    double err_p = 0.0;

                    if (nt != '-') {
                        seq_pos++;
                        err_p = phred_err(reads[i].quality[seq_pos]);
                    }

                    if (seq_pos >= 0 && seq_pos < reads[i].quality.size()) {
                        char cnt = consensus_nt[k];
                        auto consensus_info = nt_info[k][cnt];

                        double occ_ratio = double(consensus_info.occ) / double(consensus_info.total_occ);

                        // consensus is gap
                        if (cnt == '-') {
                            if (nt == '-') {
                                // gap 2 gap
                                // res_read += cnt;
                            } else {
                                // nt 2 gap (delete possible insertion)
                                if (occ_ratio >= gap_occ) {
                                    // res_read += cnt;
                                    n2g++;
                                } else {
                                    res_read += nt;
                                    res_qt += reads[i].quality[seq_pos];
                                }
                            }
                        } else {
                            if (nt == '-') {
                                // gap 2 nt (fix possible deletion)
                                if (occ_ratio >= gap_occ) {
                                    res_read += cnt;
                                    res_qt += phred_symbol(consensus_info.err);
                                    g2n++;
                                // } else {
                                    // res_read += nt;
                                }
                            } else {
                                if (nt == cnt) {
                                    // same base
                                    res_read += nt;
                                    res_qt += reads[i].quality[seq_pos];
                                } else {
                                    // sub
                                    if (occ_ratio >= min_occ && err_ratio * err_p > consensus_info.err) { // strict > to avoid subs in re-alignments
                                        res_read += cnt;
                                        res_qt += phred_symbol(consensus_info.err);
                                        n2n++;
                                    } else {
                                        res_read += nt;
                                        res_qt += reads[i].quality[seq_pos];
                                    }
                                }
                            }
                        }

                        if (seq_pos == reads[i].quality.size() - 1) {
                            seq_pos++; //end of seq
                            // break;
                        }
                    }
                }

                // res_read.erase(std::remove(res_read.begin(), res_read.end(), '-'), res_read.end());
                // {
                //     std::lock_guard<std::mutex> lock(mu);
                
                    if (res_read.size() > 0) {
                        corrected_reads.push_back(read_t{reads[i].header, res_read, "+", res_qt});
                    } else {
                        uncorrected_reads.push_back(reads[i]);
                    }
                // }
            }  
        // }));
    }

    // for (auto &&task : tasks) {
    //     task.get();
    // }

    // print_vector(consensus_nt);
    consensus_nt.erase(std::remove(consensus_nt.begin(), consensus_nt.end(), '-'), consensus_nt.end());
    // print_vector(consensus_nt);
    std::string consensus(consensus_nt.begin(), consensus_nt.end());
    // std::cout << "C: " << consensus << std::endl;
    return corrected_pack_t{-1, consensus, corrected_reads, uncorrected_reads};
}

correction_results_t correct_reads(const cluster_set_t &clusters, read_set_t &reads, double min_occ, double gap_occ, double err_ratio, int split, int min_reads, int n_threads, bool verbose, std::vector<std::string> labels) {
    std::queue<pack_to_correct_t> pending_clusters;
    int corrected = 0;
    int total_reads = 0;
    int cid = 0;

    read_set_t uncorrected_read_set;
    read_set_t corrected_read_set;
    read_set_t consensus_set;

    for (auto &tc: clusters) {           
        // int n_files = (tc.seqs.size() + split - 1) / split; // ceil(tc.seqs.size / split)
        // Avoid out of bound
        int n_files = (tc.seqs.size() - 1) / split + 1;
        int gid = tc.main_seq.gene_id; 

        // Split clusters
        for (int nf = 0; nf < n_files; nf++) {
            // int nreads_in_cluster = (tc.seqs.size() + n_files - 1 - nf) / n_files;
            int nreads_in_cluster = (tc.seqs.size() - 1 - nf) / n_files + 1;
            auto creads = read_set_t(nreads_in_cluster);

            int i = 0;
            for (int j = nf; j < tc.seqs.size(); j += n_files) {
                auto ts = tc.seqs[j];
                if (ts.rev) {
                    reads[ts.seq_id].seq = reverse_complement(reads[ts.seq_id].seq);                
                    std::reverse(reads[ts.seq_id].quality.begin(), reads[ts.seq_id].quality.end()); 
                }

                reads[ts.seq_id].header = reads[ts.seq_id].header + ",gene_cluster_" + std::to_string(gid);
                creads[i] = reads[ts.seq_id];
                ++i;
                ++total_reads;
            }

            if (creads.size() > min_reads) {
                pending_clusters.push(pack_to_correct_t{cid, creads});
            } else {
                for (int i = 0; i < creads.size(); ++i) {
                    ++corrected;
                    uncorrected_read_set.push_back(creads[i]);
                }
            }
        }
        ++cid;
    }

    std::mutex mu;
    std::vector<std::future<void>> tasks;
    auto consensi = std::vector<read_set_t>(clusters.size());

    int nf = 0;
    for (int t = 0; t < n_threads; ++t) {
        tasks.emplace_back(std::async(std::launch::async, [&nf, t, &consensi, &corrected_read_set, &uncorrected_read_set, &consensus_set, &pending_clusters, &mu, &corrected, &total_reads, gap_occ, min_occ, n_threads, &verbose, &labels] {
            while (true) {
                pack_to_correct_t pack;
                int a;

                {
                    std::lock_guard<std::mutex> lock(mu);
                    if (pending_clusters.empty()) break;

                    pack = pending_clusters.front();
                    pending_clusters.pop();
                    a = nf++;

                    if(verbose) print_progress(corrected, total_reads);
                }

                auto creads = pack.reads;
                auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(0),
                5, -4, -8, -6);

                auto graph = spoa::createGraph();
                for (int j = 0; j < creads.size(); ++j) {
                    auto alignment = alignment_engine->align(creads[j].seq, graph);
                    graph->add_alignment(alignment, creads[j].seq);
                }
                
                std::vector<std::string> msa;
                graph->generate_multiple_sequence_alignment(msa);
                
                fix_msa_ends(creads, msa);

                auto corrected_reads_pack = correct_read_pack(creads, msa, min_occ, gap_occ, 30.0, 1);
                // auto corrected_reads_pack = correct_read_pack(creads, msa, min_occ, gap_occ, 30.0, n_threads);
                auto corrected_reads = corrected_reads_pack.reads;
                auto uncorrected_reads = corrected_reads_pack.uncorrected_reads;
                {
                    std::lock_guard<std::mutex> lock(mu);
                    for (int i = 0; i < corrected_reads.size(); ++i) {
                        corrected_read_set.push_back(corrected_reads[i]);
                    }

                    // add uncorrected reads from the pack
                    for (int i = 0; i < uncorrected_reads.size(); ++i) {
                        uncorrected_read_set.push_back(uncorrected_reads[i]);
                    }

                    corrected+=creads.size();
                }
                // create new MSA with corrected reads
                sort_read_set(corrected_reads);
                graph = spoa::createGraph();

                for (int j = 0; j < corrected_reads.size(); ++j) {
                    auto alignment = alignment_engine->align(corrected_reads[j].seq, graph);
                    graph->add_alignment(alignment, corrected_reads[j].seq);
                }

                msa = std::vector<std::string>();
                graph->generate_multiple_sequence_alignment(msa);
                fix_msa_ends(corrected_reads, msa);

                // TODO: check and compare to use 1 or n_threads
                auto cv = generate_consensus_vector(corrected_reads, msa, 1);
                // auto cv = generate_consensus_vector(corrected_reads, msa, n_threads);
                // print_vector(cv.consensus_nt);
                cv.consensus_nt.erase(std::remove(cv.consensus_nt.begin(), cv.consensus_nt.end(), '-'), cv.consensus_nt.end());
                // print_vector(cv.consensus_nt);
                std::string consensus(cv.consensus_nt.begin(), cv.consensus_nt.end());

                {
                    std::lock_guard<std::mutex> lock(mu);
                    
                    // save in consensus header the gene cluster id (if applicable)
                    // the number of reads of this cluster
                    // add file labels & labels counts to consensi header
                    std::vector<std::string> labelset;     
                    std::string gid;        
                    for(auto r: creads){
                        int index = r.header.find_first_of(",");
                        int i = r.header.find_last_of(",");
                        std::string label = r.header.substr(index + 1, i - index - 1);
                        labelset.push_back(label);
                        index = r.header.find_last_of("_");
                        gid = r.header.substr(index + 1);
                    }

                    std::string label_result;
                    for(auto label: labels){
                        label_result = label_result + " " + label + ":" + std::to_string(count(labelset.begin(), labelset.end(), label));
                    }
                    consensi[pack.original_cluster_id].push_back(read_t{gid + "," + std::to_string(creads.size()) + "," + label_result, consensus, "+", std::string(consensus.size(), 'K')});

                    // // sort corrected cluster
                    // std::stable_sort(corrected_reads.begin(), corrected_reads.end(), [](read_t a, read_t b) {
                    //     return a.seq.size() > b.seq.size();
                    // });
                }
            }
        }));
    }

    for (auto &&task : tasks) {
        task.get();
    }
    
    if(verbose) print_progress(corrected, total_reads);
    std::cerr << std::endl;

    std::cerr << "Generating consensi..." << std::endl;
    cid = 0;
    for (auto& it: consensi) {
        int total_reads = 0;
        std::vector<int> label_counts(labels.size());
        std::string labels_result="";
        int gid;

        for (const auto& rit: it) {
            auto num = splitString(rit.header, ',');
            gid = std::stoi(num[0]);
            total_reads += std::stoi(num[1]);

            int index, i = 0;
            for(auto label: labels){
                if(rit.header.find(label) != std::string::npos){
                    index = rit.header.find(label);
                    std::string sub = rit.header.substr(index + 1);
                    index = sub.find_first_of(":");
                    // std::cout << label << "  " << sub.substr(index + 1) << "  " << std::stoi(sub.substr(index + 1)) << std::endl;
                    label_counts[i] += std::stoi(sub.substr(index + 1));
                }
                ++i;
            }
        }

        int  i = 0;
        for(auto label: labels){
            labels_result += label + ":" + std::to_string(label_counts[i]) + ",";
            ++i;
        }

        if (it.size() > 1) {
            auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(0),
                5, -4, -8, -6);

            auto graph = spoa::createGraph();

            for (int j = 0; j < it.size(); ++j) {
                auto alignment = alignment_engine->align(it[j].seq, graph);
                graph->add_alignment(alignment, it[j].seq);
            }
            
            //std::string consensus = graph->generate_consensus();
            std::vector<std::string> msa;
            graph->generate_multiple_sequence_alignment(msa);
            fix_msa_ends(it, msa);

            auto cv = generate_consensus_vector(it, msa, n_threads);
            cv.consensus_nt.erase(std::remove(cv.consensus_nt.begin(), cv.consensus_nt.end(), '-'), cv.consensus_nt.end());
            std::string consensus(cv.consensus_nt.begin(), cv.consensus_nt.end());
            
            if(gid != -1){
                consensus_set.push_back(read_t{"@cluster_" + std::to_string(cid) + " gene_cluster_" + std::to_string(gid) + " reads=" + std::to_string(total_reads) + " labels=" + labels_result, consensus, "+", std::string(consensus.size(), 'K')});
            } else {
                consensus_set.push_back(read_t{"@cluster_" + std::to_string(cid) + " reads=" + std::to_string(total_reads) + " labels=" + labels_result, consensus, "+", std::string(consensus.size(), 'K')});
            }
        } else {
            if (it.size() > 0) {
                if(gid != -1){
                    consensus_set.push_back(read_t{"@cluster_" + std::to_string(cid) + " gene_cluster_" + std::to_string(gid)+ " reads=" + std::to_string(total_reads) + " labels=" + labels_result, it[0].seq, "+", it[0].quality});
                }else {
                    consensus_set.push_back(read_t{"@cluster_" + std::to_string(cid) + " reads=" + std::to_string(total_reads) + " labels=" + labels_result, it[0].seq, "+", it[0].quality});
                }
            }
        }

        ++cid;
        if(verbose) print_progress(cid, consensi.size());
    }

    return correction_results_t{
        corrected_read_set,
        uncorrected_read_set,
        consensus_set
    };
}