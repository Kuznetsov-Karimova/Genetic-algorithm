#pragma once

#include <algorithm>
#include <iterator>
#include <map>
#include <optional>
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <random>
#include <numeric>
#include <iostream>


struct City {
    uint64_t number; // номер города

    double x = 0; // для координатных
    double y = 0;

    City() = default;

    explicit City(uint64_t num) {
        number = num;
    }

    City(uint64_t num, double x_coor, double y_coor) {
        number = num;
        x = x_coor;
        y = y_coor;
    }

    auto operator==(const City& other) const -> bool {
        return this->number == other.number;
    }
};


auto coord_distance(const City& a, const City& b) -> double {
    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
}


enum class MutationType { SWAP, INVERSION, DISPLACEMENT };


class Chromosome {
public:
    explicit Chromosome(const std::vector<City>& path) {
        m_full_path.assign(path.begin(), path.end());
    }

    Chromosome() = default;
    
    // mutation
    void mutation(MutationType type, double mutation_probability) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        if (dis(gen) > mutation_probability) {
            return; // Probability check
        }

        if (m_full_path.size() <= 1) {
            return;
        }

        switch (type) {
        case MutationType::SWAP: {
            std::uniform_int_distribution<int> dist(0, static_cast<int>(m_full_path.size() - 1));
            int idx1 = dist(gen);
            int idx2 = dist(gen);
            std::swap(m_full_path[idx1], m_full_path[idx2]);

            break;
        }
        case MutationType::INVERSION: {
            std::uniform_int_distribution<int> dist(0, static_cast<int>(m_full_path.size() - 1));
            int idx1 = dist(gen);
            int idx2 = dist(gen);

            if (idx1 > idx2) { std::swap(idx1, idx2); }

            std::reverse(m_full_path.begin() + idx1, m_full_path.begin() + idx2 + 1);

            break;
        }
        case MutationType::DISPLACEMENT: {
            std::uniform_int_distribution<int> dist(0, static_cast<int>(m_full_path.size() - 1));
            int idx1 = dist(gen);
            int idx2 = dist(gen);

            if (idx1 > idx2) { std::swap(idx1, idx2); }

            std::vector<City> sub_path(m_full_path.begin() + idx1, m_full_path.begin() + idx2 + 1);
            m_full_path.erase(m_full_path.begin() + idx1, m_full_path.begin() + idx2 + 1);
            int insert_pos = static_cast<int>(dist(gen) % m_full_path.size());
            m_full_path.insert(m_full_path.begin() + insert_pos, sub_path.begin(), sub_path.end());

            break;
        }
        }
    }

    [[nodiscard]] auto fitness() const -> double {
        double full_sum = 0;
        for (int i = 0; i < m_full_path.size() - 1; ++i) {
            full_sum += coord_distance(m_full_path[i], m_full_path[i+1]);
        }
        full_sum += coord_distance(m_full_path[m_full_path.size() - 1], m_full_path[0]);

        return full_sum;
    }

    [[nodiscard]] auto fitness(const std::vector<std::vector<int>>& adjacency_matrix) const -> double{
        double full_sum = 0;
        for (int i = 0; i < m_full_path.size() - 1; ++i) {
            full_sum += adjacency_matrix[m_full_path[i].number][m_full_path[i+1].number];
        }
        full_sum += adjacency_matrix[m_full_path[m_full_path.size() - 1].number][m_full_path[0].number];

        return full_sum;
    }

    std::vector<City> m_full_path; 
};


struct Parameters {
    Parameters(uint16_t gen_count, uint16_t pop_size, double prob) : generations_count(gen_count),
                                        population_size(pop_size), mutation_probability(prob) {}

    uint16_t generations_count;
    uint16_t population_size;
    double mutation_probability;
};


enum class SelectionType { ROULETTE_WHEEL, TOURNAMENT, RANK };
enum class CrossoverType { PARTIALLY_MAPPED, ORDERED, CYCLE };


class TSPPopulation {
public:
    explicit TSPPopulation(const Parameters& params, MutationType m_type, SelectionType s_type,
                        CrossoverType c_type)
        : generations_count(params.generations_count), population_size(params.population_size),
          mutation_probability(params.mutation_probability), mutation_type(m_type),
          selection_type(s_type), crossover_type(c_type)  {}

    void set_adjacency_mat(std::vector<std::vector<int>>& adjacency_matrix_info) {
        adjacency_mat = adjacency_matrix_info;
        COORD_TYPE = false;
    }

    void initialize_first_population(const std::vector<City>& cities) {
        std::random_device rd;
        std::mt19937 g(rd());

        for (uint16_t i = 0; i < population_size; ++i) {
            std::vector<City> shuffled_cities(cities);
            std::shuffle(shuffled_cities.begin(), shuffled_cities.end(), g);

            Chromosome newChromosome(shuffled_cities);
            m_population.push_back(newChromosome);
        }
    }

    void replace_population(const std::vector<Chromosome>& new_population) {
        m_population.clear();

        m_population.insert(m_population.end(), new_population.begin(), new_population.end());
    }

    // selection - return parents
    auto selection(SelectionType type) -> std::pair<Chromosome, Chromosome> {
        std::random_device rd;
        std::mt19937 gen(rd());

        switch (type) {
        case SelectionType::ROULETTE_WHEEL: {
            std::vector<double> fitness_values;
            fitness_values.reserve(m_population.size());
            for (const auto &chromo : m_population) {
                if (COORD_TYPE) { fitness_values.push_back(1.0 / chromo.fitness()); }
                else { fitness_values.push_back(1.0 / chromo.fitness(adjacency_mat)); }
            }
            std::partial_sum(fitness_values.begin(), fitness_values.end(), fitness_values.begin());

            std::uniform_real_distribution<> dis(0.0, fitness_values.back());

            auto select_one = [&]() -> Chromosome {
                double r = dis(gen);
                auto it = std::lower_bound(fitness_values.begin(), fitness_values.end(), r);
                return m_population[it - fitness_values.begin()];
            };

            return std::make_pair(select_one(), select_one());

            break;
        }

        case SelectionType::TOURNAMENT: {
            std::uniform_int_distribution<int> dist(0, population_size - 1);

            auto select_one = [&]() -> Chromosome {
                Chromosome best = m_population[dist(gen)];
                for (int i = 1; i < 3; ++i) {
                    Chromosome temp = m_population[dist(gen)];
                    if (COORD_TYPE) {
                        if (temp.fitness() < best.fitness()) {
                            best = temp;
                        }
                    }
                    else {
                        if (temp.fitness(adjacency_mat) < best.fitness(adjacency_mat)) {
                            best = temp;
                        }
                    }
                }
                return best;
            };

            return std::make_pair(select_one(), select_one());

            break;
        }

        case SelectionType::RANK: {
            std::sort(m_population.begin(), m_population.end(), [&](const Chromosome &a, const Chromosome &b) {
                if (COORD_TYPE) {return a.fitness() < b.fitness(); }
                return a.fitness(adjacency_mat) < b.fitness(adjacency_mat); });

            std::vector<int> ranks(population_size);
            std::iota(ranks.begin(), ranks.end(), 1);

            double total_rank = std::accumulate(ranks.begin(), ranks.end(), 0.0);
            std::vector<double> selection_probs(population_size);

            for (int i = 0; i < population_size; ++i) {
                selection_probs[i] = ranks[i] / total_rank;
            }

            std::partial_sum(selection_probs.begin(), selection_probs.end(), selection_probs.begin());
            std::uniform_real_distribution<> dis(0.0, selection_probs.back());

            auto select_one = [&]() -> Chromosome {
                double r = dis(gen);
                auto it = std::lower_bound(selection_probs.begin(), selection_probs.end(), r);
                return m_population[it - selection_probs.begin()];
            };

            return std::make_pair(select_one(), select_one());

            break;
        }

        return {m_population[0], m_population[1]};
    }
    }

    auto find_city_in_other_parent(int index, const std::vector<City>& parent1, const std::vector<City>& parent2) -> int {
        for (int i = 0; i < parent1.size(); ++i) {
            if (parent1[index].number == parent2[i].number) {
                return i;
            }
        }
        return -1;
    }

    // crossover - children from parents
    void crossover(Chromosome& parent1, Chromosome& parent2, CrossoverType type) {
        std::random_device rd;
        std::mt19937 gen(rd());

        switch(type) {
        case CrossoverType::PARTIALLY_MAPPED: {
            std::uniform_int_distribution<int> dist(0, static_cast<int>(parent1.m_full_path.size() - 1));
            int start = dist(gen);
            int end = dist(gen);

            if (start > end) { std::swap(start, end); }

            std::vector<City> child1(parent1.m_full_path.size());
            std::vector<City> child2(parent1.m_full_path.size());

            std::copy(parent1.m_full_path.begin() + start, parent1.m_full_path.begin() + end + 1, child1.begin() + start);
            std::copy(parent2.m_full_path.begin() + start, parent2.m_full_path.begin() + end + 1, child2.begin() + start);

            for (int i = 0; i < parent1.m_full_path.size(); ++i) {
                if (i < start || i > end) {
                    int p1_index = find_city_in_other_parent(i, parent1.m_full_path, parent2.m_full_path);
                    int p2_index = find_city_in_other_parent(i, parent2.m_full_path, parent1.m_full_path);

                    if (p1_index != -1) {
                        child1[i] = parent2.m_full_path[p1_index];
                    } else {
                        child1[i] = parent1.m_full_path[i];
                    }

                    if (p2_index != -1) {
                        child2[i] = parent1.m_full_path[p2_index];
                    } else {
                        child2[i] = parent2.m_full_path[i];
                    }
                }
            }

            parent1.m_full_path = std::move(child1);
            parent2.m_full_path = std::move(child2);

            break;
        }
        case CrossoverType::ORDERED: {
            std::uniform_int_distribution<int> dist(0, static_cast<int>(parent1.m_full_path.size() - 1));
            int start = dist(gen);
            int end = dist(gen);

            if (start > end) { std::swap(start, end); }

            auto child1 = parent1.m_full_path;
            auto child2 = parent2.m_full_path;

            child1.erase(child1.begin() + start, child1.begin() + end + 1);
            child2.erase(child2.begin() + start, child2.begin() + end + 1);

            std::vector<City> temp1;
            std::vector<City> temp2;

            std::copy_if(parent2.m_full_path.begin(), parent2.m_full_path.end(), std::back_inserter(temp1), [&child1](City& city) {
                return std::find(child1.begin(), child1.end(), city) == child1.end();
            });

            std::copy_if(parent1.m_full_path.begin(), parent1.m_full_path.end(), std::back_inserter(temp2), [&child2](City& city) {
                return std::find(child2.begin(), child2.end(), city) == child2.end();
            });

            auto pos1 = child1.begin() + start;
            child1.insert(pos1, parent1.m_full_path.begin() + start, parent1.m_full_path.begin() + end + 1);
            std::rotate(child1.begin(), child1.end() - (end + 1 - start), child1.end());

            auto pos2 = child2.begin() + start;
            child2.insert(pos2, parent2.m_full_path.begin() + start, parent2.m_full_path.begin() + end + 1);
            std::rotate(child2.begin(), child2.end() - (end + 1 - start), child2.end());

            parent1.m_full_path = std::move(child1);
            parent2.m_full_path = std::move(child2);

            break;

        }
        case CrossoverType::CYCLE: {
            auto cycle = [&](Chromosome& parent1, Chromosome& parent2) {
                std::vector<City> child(parent1.m_full_path.size());
                std::vector<bool> visited(parent1.m_full_path.size(), false);

                int index = 0;
                while (index < parent1.m_full_path.size() && !visited[index]) {
                    visited[index] = true;
                    child[index] = parent1.m_full_path[index];
                    index = static_cast<int>(std::distance(parent1.m_full_path.begin(),
                        std::find(parent1.m_full_path.begin(), parent1.m_full_path.end(), parent2.m_full_path[index])));
                }

                for (size_t i = 0; i < child.size(); ++i) {
                    if (!visited[i]) {
                        child[i] = parent2.m_full_path[i];
                    }
                }

                return child;
            };

            parent1.m_full_path = cycle(parent1, parent2);
            parent2.m_full_path = cycle(parent2, parent1);

            break;
        }
        }
    }


    void genetic_algorithm(){
        for (uint16_t generation = 0; generation < generations_count; ++generation) {
            std::vector<Chromosome> new_population;

            for (uint16_t i = 0; i < population_size / 2; ++i) {
                // selection
                auto [parent1, parent2] = selection(selection_type);

                // crossover
                Chromosome child1(parent1.m_full_path);
                Chromosome child2(parent2.m_full_path);
                crossover(child1, child2, crossover_type);

                // mutation
                child1.mutation(mutation_type, mutation_probability);
                child2.mutation(mutation_type, mutation_probability);

                new_population.push_back(std::move(child1));
                new_population.push_back(std::move(child2));
            }

            replace_population(new_population);
        }
    }

    auto get_answer() -> double {
        if (m_population.empty()) {
            std::cout << "No solution found" << std::endl;
            return 0.0;
        }

        Chromosome min_fit_chromosome = *std::min_element(this->m_population.begin(), this->m_population.end(),
                                                    [this](const Chromosome& a, const Chromosome& b) {
                                                        if (COORD_TYPE) { return a.fitness() < b.fitness(); }
                                                         return a.fitness(adjacency_mat) < b.fitness(adjacency_mat);
                                                    });

        double shortest_path = COORD_TYPE ? min_fit_chromosome.fitness() : min_fit_chromosome.fitness(adjacency_mat);
        std::cout << "PATH:" << std::endl;
        auto full_path = min_fit_chromosome.m_full_path;
        for (auto it = full_path.begin(); it != full_path.end(); ++it) {
            std::cout << it->number << " -> ";
        }
        std::cout << full_path[0].number << std::endl; // print the last city

        return shortest_path;
    }


protected:
    std::vector<Chromosome> m_population;

    bool COORD_TYPE = true;

    std::vector<std::vector<int>> adjacency_mat;

    uint16_t generations_count;
    uint16_t population_size;
    double mutation_probability;

    MutationType mutation_type;
    SelectionType selection_type;
    CrossoverType crossover_type;
};