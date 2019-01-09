#ifndef _rattle_similarity_h
#define _rattle_similarity_h

#include <vector>
#include "kmer.hpp"

struct similarity_res_t {
    std::vector<std::pair<int,int>> lis;
    int llis;
    int bases;
    std::vector<int> distances;
};

similarity_res_t calc_similarity(std::vector<kmer_match_t> common, int kmer_size);

#endif