#include "Knapsack.hpp"
#include <random>
#include <chrono>

using chromosome = std::vector<bool>;
using vector = std::vector<int>;
using Population = std::vector<chromosome*>;
using Fitness = std::vector<int>;

size_t Nampack::algorithm(bool debug) {
    return 0;
}

void nampack_app(std::vector<std::chrono::duration<double>>& all_results, bool isPrint, bool debug = false);

class Evolution {

public:

    Evolution(const std::string& path_to_values, const std::string& path_to_sizes,
    const std::string& path_to_nap_size, const std::string& path_to_res,
    int max_epoh, size_t population_size, int mutation_frequency):
    cross_rd(), cross_generator(cross_rd),
    mut_rd(), mut_generator(mut_rd) {
        
        this->max_epoh = max_epoh;
        this->population_size = population_size;
        this->mutation_frequency = mutation_frequency;


        read_file(path_to_values, prices);
        read_file(path_to_sizes, weights);
        read_file(path_to_res, m_res);
        read_file(path_to_nap_size, nap_size);

        chromosome_len = prices.size();
        best_chromosome = chromosome(chromosome_len, false);
        population = new Population(population_size, nullptr);

        gen_first_population(*this->population, population_size, weights);
        update_fitness();

        cross_dist = std::uniform_int_distribution<>(0, chromosome_len - 1);
        mut_dist = std::uniform_int_distribution<>(0, mutation_frequency - 1);
    }

    void algorithm() {

        for (int epoh_number = 1; epoh_number < max_epoh; ++epoh_number) {

            update_population();
            update_fitness();

        }

        best_chromosome = (*best_solution());
    }

    void print_chromosome(const chromosome& ch) const {

        for (int elem: ch) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;

    }

    void print_values() const {

        for (int elem: prices) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;

    }

    void print_weights() const {

        for (int elem: weights) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;

    }

    [[nodiscard]] int get_count_of_opers() const { return m_count_of_opers; } 
    [[nodiscard]] int get_nap_size() const { return nap_size; } 
    [[nodiscard]] int get_res() const { return result; }


private:

    void gen_first_population(Population population, int population_size, vector& weights) {
        
        std::random_device rd;   
        std::mt19937 generator(rd());
        std::uniform_int_distribution<> dist(0, 1);

        for(size_t chrom = 0; chrom < population_size; ++chrom) {

            population[chrom] = new chromosome(chromosome_len, false);

            for(size_t gen = 0; gen < chromosome_len; ++gen) {
                (*population[chrom]).push_back(dist(generator));
            }

        }
    }

    void update_population() {

        Population* new_population = new Population(population_size, nullptr);
        (*new_population)[0] = best_solution();

        for(size_t population_ind = 1; population_ind < population_size; ++population_ind) {

            std::pair<size_t, size_t> chroms;
            select(chroms);
            chromosome* new_chrom = crossover((*population)[chroms.first],
                                            (*population)[chroms.second]);
            mutation(*new_chrom);
            (*new_population)[population_ind] = new_chrom;
        }

        population = new_population;

    }

    void update_fitness() {
        
        Fitness* new_fitness = new Fitness(population_size, 0);

        for(size_t chrom_number = 0; chrom_number < population_size; ++chrom_number) {
            int total_weight = nap_size + 1;
            while(total_weight > nap_size) {
                int contribution = 0;
                total_weight = 0;
                vector indexes = {};
                size_t ind_size = 0;

                for(size_t gen = 0; gen < chromosome_len; ++gen) {
                    
                    if ((*(*population)[chrom_number])[gen] == true) {
                        indexes.push_back(gen);
                        ++ind_size;
                        total_weight += weights[gen];
                        contribution += prices[gen];
                    }
                }

                if (total_weight > nap_size) {

                    std::random_device rd;   
                    std::mt19937 generator(rd());
                    std::uniform_int_distribution<> dist(0, ind_size - 1);

                    int gen = dist(generator);
                    (*(*population)[chrom_number])[gen] == false;

                }

            }

            (*new_fitness)[chrom_number] = total_weight;
        
        }

        fitness = new_fitness;

    }

    chromosome* crossover(chromosome* ch_1, chromosome* ch_2) {

        chromosome* new_ch = new chromosome(chromosome_len, false);
        size_t cross_point = cross_dist(cross_generator);

        for(size_t gen = 0; gen < cross_point; ++gen) {
            (*new_ch)[gen] = (*ch_1)[gen];
        }

        for(size_t gen = cross_point; gen < chromosome_len; ++gen) {
            (*new_ch)[gen] = (*ch_2)[gen];
        }

    }

    void mutation(chromosome& ch) {

        for(size_t gen = 0; gen < chromosome_len; ++gen) {
            
            if ((mut_dist(mut_generator)) == 0) {

                if (ch[gen] == false) {
                    ch[gen] = true;
                } else {
                    ch[gen] = false;
                }
            }
        }
    }

    chromosome* best_solution() {
        
        size_t sol_index = 0;

        for(size_t chrom = 0; chrom < population_size; ++chrom) {

            if(fitness[chrom] > fitness[sol_index]) {
                sol_index = chrom;
            }

        }

        return (*population)[sol_index];
    
    }

    void select(std::pair<size_t, size_t>& res) {

    }


    int max_epoh;
    size_t population_size;
    int mutation_frequency;


    vector prices;
    vector weights;
    size_t m_res;
    size_t nap_size;

    size_t chromosome_len;
    chromosome best_chromosome = {};
    int result = 0;

    Population* population;
    Fitness* fitness;

    std::random_device cross_rd; 
    std::mt19937 cross_generator;
    std::uniform_int_distribution<> cross_dist;

    std::random_device mut_rd; 
    std::mt19937 mut_generator;
    std::uniform_int_distribution<> mut_dist;

    int m_count_of_opers = 0;

};



auto main() -> int {

    std::vector<std::chrono::duration<double>> all_results(9, static_cast<std::chrono::duration<double>>(0));

    std::vector<std::chrono::duration<double>> one_run(9);


    nampack_app(one_run, true);

    for (int j = 0; j < 9; ++j) {
        all_results[j] += one_run[j];
    }

    return 0;
}

void nampack_app(std::vector<std::chrono::duration<double>>& all_results, bool isPrint, bool debug) {


    // ALL TESTS
    std::chrono::duration<double> all_elapsed{};

    int num = 0;

    for (auto test: test_files) {
        if (isPrint) {
            std::cout << "Test: " << num << std::endl;
        }

        Evolution evol(test[0], test[1], test[2], test[3], 100, 20, 10);
        
        std::cout << "Capacity: " << evol.get_nap_size() << std::endl;
        std::cout << "Prices: ";
        evol.print_values();
        std::cout << "Weights: ";
        evol.print_weights();

        // ALGHORITM 1 TEST TIME START
        auto start = std::chrono::high_resolution_clock::now();

        evol.algorithm();

        auto end = std::chrono::high_resolution_clock::now();

        // ALGHORITM 1 TEST TIME END
        std::chrono::duration<double> elapsed = end - start;

        all_elapsed += elapsed;

        if (isPrint) {
            std::cout << "Result: " << evol.get_res() << std::endl;
            std::cout << "Time: " << elapsed.count() << " s" << std::endl;
            std::cout << "Count of opers: " << evol.get_count_of_opers() << std::endl;

            std::cout << std::endl;
        }

        all_results[num] = elapsed;

        ++num;
    }

    auto average_time = all_elapsed/(test_files.size());

    if (isPrint) {
        std::cout << "Average time: " << average_time.count() << std::endl;
    }

    all_results[num] = average_time;
}