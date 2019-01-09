#include "fasta.hpp"
#include <fstream>
#include <algorithm>
#include <iostream>

read_set_t read_fasta_file(std::string file) {
    read_set_t result;

    std::ifstream infile(file);
    std::string line;
    std::string header;

    while (std::getline(infile, line)) {
        if (line.size() == 0) continue;

        if (line[0] == '>') {
            header = line;
        } else {
            read_t r{header, line, "", ""};
            result.push_back(r);
        }
    }

    return result;
}

read_set_t read_fastq_file(std::string file) {
    read_set_t result;

    std::ifstream infile(file);
    std::string line;
    std::string header;
    std::string seq;
    std::string ann;
    std::string qt;
    int lineID = 0;

    while (std::getline(infile, line)) {
        if (lineID == 0) {
            header = line;
            ++lineID;
        } else if (lineID == 1) {
            seq = line;
            ++lineID;
        } else if (lineID == 2) {
            ann = line;
            ++lineID;
        } else if (lineID == 3) {
            qt = line;
            lineID = 0;

            read_t r{header, seq, ann, qt};
            result.push_back(r);
        }
    }

    return result;
}

bool _comp_read_set_desc(read_t a, read_t b) {
    return a.seq.size() > b.seq.size();
}

void sort_read_set(read_set_t &rs) {
    std::stable_sort(rs.begin(), rs.end(), _comp_read_set_desc);
}