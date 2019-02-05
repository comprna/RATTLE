#include <vector>
#include <future>
#include <mutex>
#include <string>
#include <algorithm>

#include "correct.hpp"
#include "utils.hpp"

read_set_t correct_reads(const read_set_t &reads, const read_set_t &aln, double min_occ, double err_ratio, int n_threads) {
    // generate consensus vector
    auto nt_info = std::vector<map_nt_info_t>(aln[0].seq.size());
    for (int i = 0; i < nt_info.size(); i++) {
        nt_info[i] = map_nt_info_t();
        nt_info[i]['A'] = pos_info_t{'A', 0.0, 0, 0};
        nt_info[i]['C'] = pos_info_t{'C', 0.0, 0, 0};
        nt_info[i]['T'] = pos_info_t{'T', 0.0, 0, 0};
        nt_info[i]['G'] = pos_info_t{'G', 0.0, 0, 0};
        nt_info[i]['-'] = pos_info_t{'-', 0.0, 0, 0};
    }

    std::mutex mu;
    
    std::vector<std::future<void>> tasks;
    for (int t = 0; t < n_threads; ++t) {
        tasks.emplace_back(std::async(std::launch::async, [t, &reads, &aln, n_threads, &mu, &nt_info] {
            std::vector<map_nt_info_t> local_nt_info = std::vector<map_nt_info_t>(aln[0].seq.size());
            for (int i = 0; i < local_nt_info.size(); i++) {
                local_nt_info[i] = map_nt_info_t();
                local_nt_info[i]['A'] = pos_info_t{'A', 0.0, 0, 0};
                local_nt_info[i]['C'] = pos_info_t{'C', 0.0, 0, 0};
                local_nt_info[i]['T'] = pos_info_t{'T', 0.0, 0, 0};
                local_nt_info[i]['G'] = pos_info_t{'G', 0.0, 0, 0};
                local_nt_info[i]['-'] = pos_info_t{'-', 0.0, 0, 0};
            }

            for (int i = t; i < reads.size(); i+=n_threads) {
                auto read_aln = aln[i];
                int seq_pos = -1;

                for (int k = 0; k < read_aln.seq.size(); k++) {
                    char nt = read_aln.seq[k];
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
                        }
                    }
                }
            }

            std::lock_guard<std::mutex> lock(mu);
            for (int k = 0; k < aln[0].seq.size(); ++k) {
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
    auto consensus_nt = std::vector<char>(aln[0].seq.size());
    for (int k = 0; k < aln[0].seq.size(); ++k) {
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

        consensus_nt[k] = max_nt;
    }

    // correct aln
    auto corrected_reads = read_set_t(reads.size());

    tasks.clear();
    for (int t = 0; t < n_threads; ++t) {
        tasks.emplace_back(std::async(std::launch::async, [t, &reads, &aln, n_threads, &mu, &nt_info, &consensus_nt, min_occ, err_ratio, &corrected_reads] {
            for (int i = t; i < reads.size(); i+=n_threads) {
                auto read_aln = aln[i];
                int seq_pos = -1;
                std::string res_read;
                std::string res_qt;

                int n2g = 0;
                int g2n = 0;
                int n2n = 0;

                for (int k = 0; k < read_aln.seq.size(); ++k) {
                    char nt = read_aln.seq[k];
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
                                res_read += cnt;
                            } else {
                                // nt 2 gap (delete possible insertion)
                                if (occ_ratio >= min_occ) {
                                    res_read += cnt;
                                    n2g++;
                                } else {
                                    res_read += nt;
                                    res_qt += reads[i].quality[seq_pos];
                                }
                            }
                        } else {
                            if (nt == '-') {
                                if (occ_ratio >= min_occ) {
                                    res_read += cnt;
                                    res_qt += phred_symbol(consensus_info.err);
                                    g2n++;
                                } else {
                                    res_read += nt;
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
                        }
                    }
                }

                res_read.erase(std::remove(res_read.begin(), res_read.end(), '-'), res_read.end());
                corrected_reads[i] = read_t{reads[i].header, res_read, "+", res_qt};
            }  
        }));
    }

    for (auto &&task : tasks) {
        task.get();
    }

    return corrected_reads;
}