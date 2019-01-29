#include "cluster.hpp"
#include "kmer.hpp"
#include "similarity.hpp"
#include "utils.hpp"

#include <future>
#include <mutex>
#include <string>
#include <iostream>

cseq_t cluster_together(const read_set_t &reads, const std::vector<std::vector<kmer_t>> &kmers, const std::vector<std::vector<kmer_t>> &rev_kmers, const std::vector<kmer_bv_t> &bv_kmers, const std::vector<kmer_bv_t> &rev_bv_kmers, int i, int j, int kmer_size, double t_s, double t_v, double bv_threshold) {
    // filter by kmer bv intersection score
    auto bv_common = (bv_kmers[i] & bv_kmers[j]).count();
    auto rev_bv_common = (bv_kmers[i] & rev_bv_kmers[j]).count();

    double mmax = std::max(bv_kmers[i].count(), bv_kmers[j].count());
    double bv_score = std::max(bv_common/mmax, rev_bv_common/mmax);

    if (bv_common/mmax >= bv_threshold) {
        auto common = get_common_kmers(kmers[i], kmers[j]);
        auto sim = calc_similarity(common, kmer_size);

        // normalize scores
        double mn = std::min(reads[i].seq.size(), reads[j].seq.size());
        double norm_score = double(sim.bases)/mn;

        if (norm_score >= t_s) {
            // normal strand
            if (var(sim.distances) < t_v) {
                return cseq_t{j, false};
            }
        }
    }
    
    if (rev_bv_common/mmax >= bv_threshold) {
        auto rev_common = get_common_kmers(kmers[i], rev_kmers[j]);
        auto rev_sim = calc_similarity(rev_common, kmer_size);

        // normalize scores
        double mn = std::min(reads[i].seq.size(), reads[j].seq.size());
        double rev_norm_score = double(rev_sim.bases)/mn;

        if (rev_norm_score >= t_s) {
            // reverse strand
            if (var(rev_sim.distances) < t_v) {
                return cseq_t{j, true};
            }
        }
    }

    return cseq_t{-1, false};
}

cluster_set_t cluster_reads(const read_set_t &reads, int kmer_size, double t_s, double t_v, double bv_threshold, double min_bv_threshold, double bv_falloff, int min_reads_cluster, int n_threads) {
    cluster_set_t clusters;
    std::mutex mu;
    auto already_clustered = std::vector<bool>(reads.size(), false);

    // extract kmers from reads
    std::vector<std::vector<kmer_t>> kmers(reads.size());
    std::vector<std::vector<kmer_t>> rev_kmers(reads.size());

    std::vector<kmer_bv_t> bv_kmers(reads.size());
    std::vector<kmer_bv_t> rev_bv_kmers(reads.size());
    
    std::vector<std::future<void>> tasks;
    for (int t = 0; t < n_threads; ++t) {
        tasks.emplace_back(std::async(std::launch::async, [t, &reads, n_threads, kmer_size, &kmers, &rev_kmers, &bv_kmers, &rev_bv_kmers] {
            for (int i = t; i < reads.size(); i+=n_threads) {
                read_kmers_t k1 = extract_kmers_from_read(reads[i].seq, kmer_size);

                kmers[i] = k1.list_forward;
                rev_kmers[i] = k1.list_reverse;
                bv_kmers[i] = k1.bv_forward;
                rev_bv_kmers[i] = k1.bv_reverse;
            }
        }));
    }

    for (auto &&task : tasks) {
        task.get();
    }

    // create initial clusters
    cluster_set_t initial_clusters;
    for (int i = 0; i < reads.size(); ++i) {
        std::cerr << i << std::endl;

        if (already_clustered[i]) {
            continue;
        }

        std::vector<cseq_t> cseqs;
        
        cseqs.push_back(cseq_t{i, false});
        already_clustered[i] = true;

        for (int t = 0; t < n_threads; ++t) {
            tasks.emplace_back(std::async(std::launch::async, [t, i, t_v, t_s, &mu, &kmers, &rev_kmers, &bv_kmers, &rev_bv_kmers, &reads, n_threads, kmer_size, &already_clustered, &cseqs, bv_threshold] {
                for (int j = i+1+t; j < reads.size(); j+=n_threads) {
                    if (already_clustered[j]) {
                        continue;
                    }

                    // CHECK TODO: Try saving cseqs to a local thread vector and then lock at the end
                    auto sinfo = cluster_together(reads, kmers, rev_kmers, bv_kmers, rev_bv_kmers, i, j, kmer_size, t_s, t_v, bv_threshold);
                    if (sinfo.seq_id != -1) {
                        std::lock_guard<std::mutex> lock(mu);
                        already_clustered[sinfo.seq_id] = true;
                        cseqs.push_back(sinfo);
                    }
                }
            }));
        }

        for (auto &&task : tasks) {
            task.get();
        }

        // create cluster
        cluster_t cluster;
        cluster.main_seq = cseq_t{i, false};
        cluster.seqs = cseqs;
        
        initial_clusters.push_back(cluster);
    }

    clusters = initial_clusters;

    // iteratively merge clusters
    double current_bv_threshold = bv_threshold - bv_falloff;
    bool last = false;
    while (current_bv_threshold >= min_bv_threshold || last) {
        cluster_set_t tmp_clusters;
        already_clustered = std::vector<bool>(clusters.size(), false);

        for (int i = 0; i < clusters.size(); ++i) {
            std::cerr << i << std::endl;

            if (already_clustered[i]) {
                continue;
            }
            already_clustered[i] = true;
            std::vector<cseq_t> clusters_to_merge;
            std::vector<std::future<void>> tasks;

            // NOTE: We are storing cluster IDs, not read IDs!!!!
            clusters_to_merge.push_back(cseq_t{i, false});
            for (int t = 0; t < n_threads; ++t) {
                tasks.emplace_back(std::async(std::launch::async, [t, i, t_v, t_s, &clusters, &mu, &kmers, &rev_kmers, &bv_kmers, &rev_bv_kmers, &reads, n_threads, kmer_size, &already_clustered, &clusters_to_merge, current_bv_threshold] {
                    for (int j = i+1+t; j < clusters.size(); j+=n_threads) {
                        if (already_clustered[j]) {
                            continue;
                        }

                        // CHECK TODO: Try saving clusters_to_merge to a local thread vector and then lock at the end
                        auto sinfo = cluster_together(reads, kmers, rev_kmers, bv_kmers, rev_bv_kmers, clusters[i].main_seq.seq_id, clusters[j].main_seq.seq_id, kmer_size, t_s, t_v, current_bv_threshold);
                        if (sinfo.seq_id != -1) {
                            std::lock_guard<std::mutex> lock(mu);
                            already_clustered[j] = true;
                            clusters_to_merge.push_back(cseq_t{j, sinfo.rev});
                        }
                    }
                }));
            }
            
            for (auto &&task : tasks) {
                task.get();
            }

            // create new merged clusters
            cluster_t cluster;

            // get main seq of new cluster
            cseq_t new_main_seq = clusters_to_merge[0];
            for (auto c : clusters_to_merge) {
                auto c_a = clusters[new_main_seq.seq_id];
                auto c_b = clusters[c.seq_id];

                if (reads[c_b.main_seq.seq_id].seq.size() > reads[c_a.main_seq.seq_id].seq.size()) {
                    new_main_seq = c;
                }
            }
            cluster.main_seq = clusters[new_main_seq.seq_id].main_seq;

            // merge all clusters
            for (auto c : clusters_to_merge) {
                // get old cluster seqs
                auto old_cluster = clusters[c.seq_id];

                for (auto s : old_cluster.seqs) {
                    if (c.rev != new_main_seq.rev) {
                        s.rev = !s.rev;
                    }

                    cluster.seqs.push_back(s);
                }
            }
            tmp_clusters.push_back(cluster);
            // if (cseqs.size() >= min_reads_cluster) {
            // }
        }

        clusters = tmp_clusters;
        std::cerr << "Iteration " << current_bv_threshold << " complete" << std::endl;
        
	    if (last) break;

        current_bv_threshold -= bv_falloff;
        if (current_bv_threshold < min_bv_threshold && !last) {
            last = true;
            current_bv_threshold = 0.0;
        }
    }

    return clusters;
}
