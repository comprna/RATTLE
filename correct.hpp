#ifndef _rattle_correct_h
#define _rattle_correct_h

#include <unordered_map>
#include <string>
#include "fasta.hpp"
#include "cluster.hpp"

struct pos_info_t {
    char nt;
    double err;
    int occ;
    int total_occ;
};

typedef std::unordered_map<char, pos_info_t> map_nt_info_t;
typedef std::vector<std::string> msa_t;

struct pack_to_correct_t {
    int original_cluster_id;
    read_set_t reads;
};

struct corrected_pack_t {
    int original_cluster_id;
    std::string consensus;
    read_set_t reads;
    read_set_t uncorrected_reads;
};

struct correction_results_t {
    read_set_t corrected;
    read_set_t uncorrected;
    read_set_t consensi;
};

struct consensus_vector_t {
    std::vector<map_nt_info_t> nt_info;
    std::vector<char> consensus_nt;
};

consensus_vector_t generate_consensus_vector(const read_set_t &reads, const msa_t &aln, int n_threads);
corrected_pack_t correct_read_pack(const read_set_t &reads, const msa_t &aln, double min_occ, double gap_occ, double err_ratio, int n_threads);
correction_results_t correct_reads(const cluster_set_t &clusters, read_set_t &reads, double min_occ, double gap_occ, double err_ratio, int split, int min_reads, int n_threads, bool verbose, std::vector<std::string> labels);
void fix_msa_ends(read_set_t &reads, msa_t &aln);
std::vector<std::string> splitString(std::string str, char delimiter);

#endif