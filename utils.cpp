#include "utils.hpp"
#include <iostream>
#include <math.h>
#include <sstream>

char phred_symbol(double p) {
    return -10 * log10(p) + 33;
}

double phred_err(char c) {
    double q = c - 33;
    return pow(10.0, -q/10.0);
}

std::string reverse_complement(std::string seq) {      
    std::string res = seq;
    int len = res.length();

    for (int i = 0; i < len; ++i) {
        res[i] = base_complements.find(seq[len-1-i])->second;
    }

    return res;
}

double mean(std::vector<int> s) {
    double res = 0.0;

    for (int n : s) {
        res += double(n);
    }

    return res / double(s.size());
}

double var(std::vector<int> s) {
    // corrected two-pass algorithm (1.7), from "Algorithms for computing
	// the sample variance: Analysis and recommendations" by Chan, Tony F., Gene H. Golub,
    // and Randall J. LeVeque.

    if (s.size() == 0) return 0;
    
    double ss = 0.0;
    double compensation = 0.0;

    double m = mean(s);

    for (int n : s) {
        double d = n - m;
        ss += d*d;
        compensation += d;
    }

    return (ss - compensation*compensation/double(s.size())) / double(s.size()-1);
}

void print_progress(int a, int b) {
    double progress = double(a)/double(b);
    int width = 80;

    std::cerr << "[";
    int pos = width * progress;
    
    for (int i = 0; i < width; ++i) {
        if (i < pos) std::cerr << "=";
        else if (i == pos) std::cerr << ">";
        else std::cerr << " ";
    }

    std::cerr << "] " << a << "/" << b << " (" << progress * 100.0 << "%)";
    if (a == b) std::cerr << std::endl;
    else std::cerr << "\r";
    
    std::cerr.flush();
}

std::vector<std::string> split(const std::string &s, char sep) {
    std::vector<std::string> res;
    std::stringstream ss(s);
    std::string item;

    while (getline(ss, item, sep))
        res.push_back(item);

    return res;
}
