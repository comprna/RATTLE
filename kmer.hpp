#ifndef _rattle_kmer_h
#define _rattle_kmer_h

#include <string>
#include <stdint.h>
#include <vector>
#include <unordered_map>
#include <bitset>

typedef std::pair<uint32_t, int> kmer_t; 
typedef std::pair<int, int> kmer_match_t; 

// bitvectors
const int KMER_BV_SIZE = 6;
const int BV_SIZE = 4<<(2*(KMER_BV_SIZE-1));
typedef std::bitset<BV_SIZE> kmer_bv_t;

struct read_kmers_t {
    std::vector<kmer_t> list_forward;
    std::vector<kmer_t> list_reverse;
    kmer_bv_t bv_forward;
    kmer_bv_t bv_reverse;
};

const std::unordered_map<char, uint32_t> base_to_int = {
    {'A', 0},
    {'C', 1},
    {'T', 2},
    {'U', 2},
    {'G', 3}
};

inline uint32_t hash_kmer(std::string km) {
    uint32_t res = 0;
    for (int i = 0; i < km.length(); ++i) {
        res = res<<2 | base_to_int.find(km[i])->second;
    }

    return res;
}

read_kmers_t extract_kmers_from_read(std::string read, int kmer_size, bool both_strands);
std::vector<kmer_match_t> get_common_kmers(const std::vector<kmer_t> &k1, const std::vector<kmer_t> &k2);

#endif