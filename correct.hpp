#ifndef _rattle_correct_h
#define _rattle_correct_h

#include <unordered_map>
#include "fasta.hpp"

struct pos_info_t {
    char nt;
    double err;
    int occ;
    int total_occ;
};

struct pack_to_correct_t {
    int original_cluster_id;
    read_set_t reads;
};

struct corrected_pack_t {
    int original_cluster_id;
    std::string consensus;
    read_set_t reads;
};

typedef std::unordered_map<char, pos_info_t> map_nt_info_t;

corrected_pack_t correct_reads(const read_set_t &reads, const read_set_t &aln, double min_occ, double gap_occ, double err_ratio, int n_threads);
#endif