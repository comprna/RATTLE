#ifndef _rattle_cluster_h
#define _rattle_cluster_h

#include "fasta.hpp"
#include "hps/src/hps.h"

#include <string>
#include <vector>

struct cseq_t {
    int seq_id;
    bool rev;

    template <class B>
    void serialize(B& buf) const {
        buf << seq_id << rev;
    }

    template <class B>
    void parse(B& buf) {
        buf >> seq_id >> rev;
    }
};

struct cluster_t {
    cseq_t main_seq;
    std::vector<cseq_t> seqs;

    template <class B>
    void serialize(B& buf) const {
        buf << main_seq << seqs;
    }

    template <class B>
    void parse(B& buf) {
        buf >> main_seq >> seqs;
    }
};


typedef std::vector<cluster_t> cluster_set_t;

cluster_set_t cluster_reads(const read_set_t &reads, int kmer_size, double t_s, double t_v, double bv_threshold, double min_bv_threshold, double bv_falloff, int min_reads_cluster, bool use_hc, double repr_percentile, int n_threads);    

#endif