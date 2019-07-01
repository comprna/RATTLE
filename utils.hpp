#ifndef _rattle_utils_h
#define _rattle_utils_h

#include <string>
#include <unordered_map>
#include <vector>
#include <random>

const std::unordered_map<char, char> base_complements = {
    {'A', 'T'},
    {'C', 'G'},
    {'T', 'A'},
    {'G', 'C'}
};

const std::string str_chars = "0123456789abcdefghijklmnopqrstuvwxyz";

char phred_symbol(double p);
double phred_err(char c);
std::string reverse_complement(std::string seq);
double mean(std::vector<int> s);
double var(std::vector<int> s);
std::string random_str(std::default_random_engine eng, int sz);
void print_progress(int a, int b);

#endif