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

typedef std::unordered_map<char, pos_info_t> map_nt_info_t;

read_set_t correct_reads(const read_set_t &reads, const read_set_t &aln, double min_occ, double err_ratio, int n_threads);
#endif