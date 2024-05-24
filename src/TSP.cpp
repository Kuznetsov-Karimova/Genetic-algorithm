#include "TSP.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

void cities_from_file(std::vector<City>& cities, std::string& edge_weight_type,
                    std::vector<std::vector<int>>& adjacency_mat, const std::string& file_name) {
    std::ifstream file(file_name);

    if (!file) {
        std::cerr << "Error: can't open file";
    }

    std::string line;
    std::string token;

    int count_of_cities;

    while (getline(file, line)) {
        if (line.find("DIMENSION") != std::string::npos) {
            std::istringstream iss(line);
            iss >> token >> token >> count_of_cities;
            break;
        }
    }

    std::string edge_weight_format;

    while (getline(file, line)) {
        if (line.find("EDGE_WEIGHT_TYPE") != std::string::npos) {
            std::istringstream iss(line);
            iss >> token >> token >> edge_weight_type;

            if (edge_weight_type == "ATT" || edge_weight_type == "EUC2D" || edge_weight_type == "EUC_2D") {
                edge_weight_type = "coord";
                break;
            }
            edge_weight_type = "matrix";
        }
        if (line.find("EDGE_WEIGHT_FORMAT") != std::string::npos) {
            std::istringstream iss(line);
            iss >> token >> token >> edge_weight_format;
            break;
        }
    }

    if (edge_weight_type == "coord") {
        int index;
        double x;
        double y;
        while (getline(file, line)) {
            if (line == "NODE_COORD_SECTION") {
                while (getline(file, line) && line != "EOF") {
                    std::istringstream iss(line);

                    iss >> index >> x >> y;
                    cities.emplace_back(index, x, y);
                }
                break;
            }
        }
    } else {
        adjacency_mat.resize(count_of_cities + 1);
        for (int city_num = 1; city_num < count_of_cities + 1; ++city_num) {
            cities.emplace_back(city_num);
            adjacency_mat[city_num].resize(count_of_cities + 1);
        }
        if (edge_weight_format == "FULL_MATRIX") {
            while (getline(file, line)) {
                if (line == "EDGE_WEIGHT_SECTION") {
                    break;
                }
            }
            for (int i = 1; i < count_of_cities + 1; ++i) {
                for (int j = 1; j < count_of_cities + 1; ++j) {
                    file >> adjacency_mat[i][j];
                }
            }
        } else {
            while (getline(file, line)) {
                if (line == "EDGE_WEIGHT_SECTION") {
                    break;
                }
            }
            int distance = 0;
            for (int i = 1; i < count_of_cities + 1; ++i) {
                for (int j = 1; j <= i ; ++j) {
                    file >> distance;
                    adjacency_mat[i][j] = distance;
                    adjacency_mat[j][i] = distance;
                }
            }
        }
    }

    file.close();
}


void gen_alg_test(const std::string& file_name) {
    std::vector<City> cities;
    std::string edge_weight_type;
    std::vector<std::vector<int>> adjacency_matrix;

    const std::string TSP_dir_path = "./benchmarks/TSP/";
    std::string path_to_file = TSP_dir_path + file_name;
    cities_from_file(cities, edge_weight_type, adjacency_matrix, path_to_file);

    double prob = 0.01; // Вероятность мутации
    // std::vector<uint16_t> gen_counts = { 100, 200, 500, 1000, 2000, 10000}; // Количество поколений
    // std::vector<uint16_t> pop_sizes = {50, 100, 200, 500}; // Размер популяции

    // fisrt - gen_count, second - population_size
    std::vector<std::pair<uint16_t, uint16_t>> params = {{100,50},
                                                        {200, 50},
                                                        {500,50},
                                                        {500, 100},
                                                        {1000, 50},
                                                        {1000, 100},
                                                        {1000, 200},
                                                        {2000, 50},
                                                        {2000, 100},
                                                        {2000, 200},
                                                        {5000, 50},
                                                        {5000, 100},
                                                        {5000, 200},
                                                        {10000, 50},
                                                        {10000, 100},
                                                        {10000,200},
                                                        };

    std::tuple<uint16_t, uint16_t> best_params;
    double best_answer = 100000000;

    // std::vector<MutationType> mut_types = {MutationType::SWAP, MutationType::DISPLACEMENT, MutationType::INVERSION};
    // std::vector<CrossoverType> cros_types = {CrossoverType::PARTIALLY_MAPPED, CrossoverType::ORDERED, CrossoverType::CYCLE};
    // std::vector<SelectionType> selec_types = {SelectionType::ROULETTE_WHEEL, SelectionType::TOURNAMENT, SelectionType::RANK};
    for (const auto& param: params) {
        // for (const auto& m_type : mut_types) {
            // for (const auto& s_type: selec_types) {
                // for (const auto& c_type: cros_types) {

        std::cout << "params: " << "gen_count = " << param.first << " ; population_size = " << param.second << std::endl;
        // std::cout << "algs: " << "mutation - " << m_type << " ; selection - " << s_type << " ; crossover - " << c_type << std::endl;
        TSPPopulation tsp_population(Parameters(param.first, param.second, prob), 
        MutationType::DISPLACEMENT, SelectionType::TOURNAMENT, CrossoverType::PARTIALLY_MAPPED);

        if (edge_weight_type == "matrix") {
            tsp_population.set_adjacency_mat(adjacency_matrix);
        }

        auto start = std::chrono::high_resolution_clock::now();
        tsp_population.initialize_first_population(cities);

        tsp_population.genetic_algorithm();
        auto end = std::chrono::high_resolution_clock::now();

        double answer = tsp_population.get_answer();
        std::cout << "Shortest path = " << answer << std::endl;

        std::chrono::duration<double> elapsed = end - start;
        std::cout << "TIME: " << elapsed << std::endl;

        if (answer < best_answer) {
            best_answer = answer;
            
            std::get<0>(best_params) = param.first; 
            std::get<1>(best_params) = param.second;
        }

        std::cout << std::endl;
    // }
    // }
    // }
    }


    std::cout << std::endl;
    std::cout << "BEST ANSWER = " << best_answer << std::endl;
    std::cout << "PARAMS:" << std::endl;
    std::cout << "generations = " << std::get<0>(best_params) << std::endl;
    std::cout << "population size = " << std::get<1>(best_params) << std::endl;
}

std::vector<std::string> files = {//"a280.txt",
                                //"att48.txt",
                                "bays29.txt",
                                "ch150.txt",
                                "fl417.txt",
                                "gr17.txt"
                                };

auto main() -> int {

    for (const auto& file : files) {
        std::cout << file << std::endl;
        gen_alg_test(file);
        std::cout << std::endl;
    }

    // TSPPopulation tsp_population(Parameters(gen_count, pop_size, prob), 
    //                     m_type, s_type, c_type);

    // if (edge_weight_type == "matrix") {
    //     tsp_population.set_adjacency_mat(adjacency_matrix);
    // }

    // auto start = std::chrono::high_resolution_clock::now();
    // tsp_population.initialize_first_population(cities);

    // tsp_population.genetic_algorithm();
    // auto end = std::chrono::high_resolution_clock::now();

    // tsp_population.print_answer();

    // std::chrono::duration<double> elapsed = end - start;
    // std::cout << "TIME: " << elapsed << std::endl;

    return 0;
}