#include "similarity.hpp"
#include <iostream>

similarity_res_t calc_similarity(std::vector<kmer_match_t> common, int kmer_size) {
    auto p = std::vector<int>(common.size());
    auto m = std::vector<int>(common.size()+1);

    int l = 0;

    for (int i = 0; i < common.size(); ++i) {
        int lo = 1;
        int hi = l;

        while (lo <= hi) {
            int mid = ((lo+hi) + 1) / 2; // ceil((lo+hi)/2)

            if (common[m[mid]].second < common[i].second) {
                lo = mid+1;
            } else {
                hi = mid-1;
            }
        }

        int new_l = lo;
        p[i] = m[new_l-1];
        m[new_l] = i;

        if (new_l > l) {
            l = new_l;
        }
    }
        
    similarity_res_t res;
    res.bases = 0;

    if (l > 0) {
        auto s = std::vector<kmer_match_t>(l);
        int k = m[l];

        // get LIS
        for (int i = l - 1; i >= 0; --i) {
            s[i] = common[k];
            k = p[k]; 
        }

        // calc bases
        int bases = 0;
        int hc_bases = 0;
        auto finalLIS = std::vector<kmer_match_t>();
        auto distances = std::vector<int>();

        for (int i = 0; i < s.size(); ++i) {
            if (i > 0) {
                if (
                    (s[i].first - finalLIS[finalLIS.size()-1].first < kmer_size &&
                    s[i].second - finalLIS[finalLIS.size()-1].second < kmer_size) ||
                    (s[i].first - finalLIS[finalLIS.size()-1].first >= kmer_size &&
                    s[i].second - finalLIS[finalLIS.size()-1].second >= kmer_size)
                ) {
                    bases += kmer_size;

                    auto ex = kmer_size - (s[i].second - s[i-1].second);
                    if (ex > 0) {
                        bases -= ex;
                    }

                    finalLIS.push_back(s[i]);

                    auto fls = finalLIS.size();
                    int dist = (finalLIS[fls-1].second - finalLIS[fls-2].second)-(finalLIS[fls-1].first - finalLIS[fls-2].first);
                    distances.push_back(dist);

                    if (dist < 10) {
                        hc_bases += kmer_size;
                        if (ex > 0) {
                            hc_bases -= ex;
                        }
                    }
                }
            } else {
                finalLIS.push_back(s[i]);
                bases += kmer_size;
                hc_bases += kmer_size;
            }
        }

        // std::cout << "bases" << bases << std::endl;

        res.lis = finalLIS;
        res.llis = finalLIS.size();
        res.bases = bases;
        res.hc_bases = hc_bases;
        res.distances = distances;
    }

    return res;
}