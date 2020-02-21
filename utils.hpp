#ifndef _rattle_utils_h
#define _rattle_utils_h

#include <string>
#include <unordered_map>
#include <vector>

const std::unordered_map<char, char> base_complements = {
    {'A', 'T'},
    {'C', 'G'},
    {'T', 'A'},
    {'G', 'C'},
    {'U', 'A'}
};

const std::string str_chars = "0123456789abcdefghijklmnopqrstuvwxyz";

char phred_symbol(double p);
double phred_err(char c);
std::string reverse_complement(std::string seq);
double mean(std::vector<int> s);
double var(std::vector<int> s);
void print_progress(int a, int b);
std::vector<std::string> split(const std::string &s, char sep);

#endif