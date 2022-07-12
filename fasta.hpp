#ifndef _rattle_fasta_h
#define _rattle_fasta_h

#include <string>
#include <vector>

struct read_t {
    std::string header;
    std::string seq;
    std::string ann;
    std::string quality;
};

typedef std::vector<read_t> read_set_t;

std::string unzip_file(std::string filename, int index);
read_set_t read_fasta_file(std::string file, std::string sample_id, int index, bool raw, int lower_len, int upper_len);
read_set_t read_fastq_file(std::string file, std::string sample_id, int index, bool raw, int lower_len, int upper_len);
read_set_t read_fasta_file(std::string file, std::string sample_id);
read_set_t read_fastq_file(std::string file, std::string sample_id);
read_set_t read_fastq_file(std::string file);
void write_fastq_file(const read_set_t &reads, std::string file);
void sort_read_set(read_set_t &rs);

#endif