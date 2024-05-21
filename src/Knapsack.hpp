#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

const std::vector<std::vector<std::string>> test_files = {
    {"./benchmarks/p01_p.txt", "./benchmarks/p01_w.txt", "./benchmarks/p01_c.txt", "./benchmarks/p01_s.txt"},
    {"./benchmarks/p02_p.txt", "./benchmarks/p02_w.txt", "./benchmarks/p02_c.txt", "./benchmarks/p02_s.txt"},
    {"./benchmarks/p03_p.txt", "./benchmarks/p03_w.txt", "./benchmarks/p03_c.txt", "./benchmarks/p03_s.txt"},
    {"./benchmarks/p04_p.txt", "./benchmarks/p04_w.txt", "./benchmarks/p04_c.txt", "./benchmarks/p04_s.txt"},
    {"./benchmarks/p05_p.txt", "./benchmarks/p05_w.txt", "./benchmarks/p05_c.txt", "./benchmarks/p05_s.txt"},
    {"./benchmarks/p06_p.txt", "./benchmarks/p06_w.txt", "./benchmarks/p06_c.txt", "./benchmarks/p06_s.txt"},
    {"./benchmarks/p07_p.txt", "./benchmarks/p07_w.txt", "./benchmarks/p07_c.txt", "./benchmarks/p07_s.txt"},
    };

void read_file(const std::string& path_to_file, std::vector<int>& vec) {

    std::ifstream in;
    in.open(path_to_file);

    int temp_int;

    if (!in.is_open()) {
        std::cout << "incorrect file name" << std::endl;
    } else {

        while(!in.eof()) {
            in >> temp_int;
            if(!in.eof()) {
                vec.push_back(temp_int);
            }
        }

        in.close();
    }
}

void read_file(const std::string& path_to_file, size_t& res) {

    std::ifstream in;
    in.open(path_to_file);

    if (!in.is_open()) {
        std::cout << "incorrect file name" << std::endl;
    } else {
        in >> res;
        in.close();
    }
}

class Nampack {
public:

    Nampack(const std::string& path_to_values, const std::string& path_to_sizes,
    const std::string& path_to_nap_size, const std::string& path_to_res) {

        read_file(path_to_values, m_prices);
        read_file(path_to_sizes, m_weight);
        read_file(path_to_res, m_res);
        read_file(path_to_nap_size, m_nap_size);

        m_count_of_sub = m_prices.size();
        m_res_subs = std::vector<int>(m_count_of_sub, 0);
    }

    void print_values() const {

        for (int elem: m_prices) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;

    }

    void print_weights() const {

        for (int elem: m_weight) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;

    }

    void print_res() const {

        for (int elem: m_res_subs) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;

    }

    size_t algorithm(bool debug = false);

    [[nodiscard]] int get_count_of_opers() const { return m_count_of_opers; } 
    [[nodiscard]] int get_nap_size() const { return m_nap_size; } 

private:
    std::vector<int> m_prices;
    std::vector<int> m_weight;
    size_t m_res;
    size_t m_nap_size;

    size_t m_count_of_sub;
    std::vector<int> m_res_subs = {};
 
    int m_count_of_opers = 0;
};