#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <map>
#include <chrono>  // For timing

using namespace std;

pair<vector<vector<int>>, vector<string>> load_csv(const string& filename) {
    vector<vector<int>> data;
    vector<string> gene_list;
    ifstream file(filename);
    string line;

    // Read the header and save gene names
    if (getline(file, line)) {
        stringstream ss(line);
        string gene_name;

        while (getline(ss, gene_name, ',')) {
            gene_list.push_back(gene_name);
        }
    }

    // Read the rest of the data
    while (getline(file, line)) {
        vector<int> row;
        stringstream ss(line);
        string cell;

        while (getline(ss, cell, ',')) {
            row.push_back(stoi(cell));
        }

        data.push_back(row);
    }

    return {data, gene_list};
}

// Function: Get gene names by indices
std::vector<std::string> get_genes_by_indices(const std::vector<std::string>& gene_list, const std::vector<int>& indices) {
    std::vector<std::string> genes;

    for (int index : indices) {
        if (index >= 0 && index < gene_list.size()) {
            genes.push_back(gene_list[index]);
        } else {
            std::cerr << "Index out of bounds: " << index << std::endl;
            genes.push_back("UNKNOWN_GENE");
        }
    }

    return genes;
}


std::string get_cancer_type(const std::string& filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
}

int choose_initial_column(const vector<vector<int>>& df) {
    int P = df.size();
    int best_col = -1;
    double min_diff = numeric_limits<double>::infinity();

    for (int col = 0; col < df[0].size(); ++col) {
        int count_ones = 0;
        for (const auto& row : df) {
            count_ones += row[col];
        }
        double diff = abs(count_ones - P / 2.0);
        if (diff < min_diff) {
            min_diff = diff;
            best_col = col;
        }
    }

    return best_col;
}

int best_column(const vector<vector<int>>& data, const vector<vector<int>>& split_list, const set<int>& unselected_columns) {
    double min_diff = numeric_limits<double>::infinity();
    int best_col = -1;

    for (int column : unselected_columns) {
        double sum = 0;
        for (const auto& part : split_list) {
            int p0 = 0, p1 = 0;
            for (int patient : part) {
                if (data[patient][column] == 0) {
                    p0++;
                } else {
                    p1++;
                }
            }
            sum += abs(p0 - p1);
        }
        if (sum < min_diff) {
            min_diff = sum;
            best_col = column;
        }
    }

    return best_col;
}

vector<vector<int>> split(const vector<vector<int>>& partition, const vector<vector<int>>& data, int column) {
    vector<vector<int>> new_partition;
    for (const auto& part : partition) {
        vector<int> p0, p1;
        for (int patient : part) {
            if (data[patient][column] == 0) {
                p0.push_back(patient);
            } else {
                p1.push_back(patient);
            }
        }
        if (!p0.empty()) new_partition.push_back(p0);
        if (!p1.empty()) new_partition.push_back(p1);
    }
    return new_partition;
}

bool at_least_one(const vector<vector<int>>& partition) {
    if (partition.size() == 1 && partition[0].size() == 1) {
        return false;
    }
    for (const auto& part : partition) {
        if (part.size() > 1) {
            return true;
        }
    }
    return false;
}

vector<vector<int>> remove_small_partition(vector<vector<int>>& partition) {
    partition.erase(remove_if(partition.begin() + 1, partition.end(),
                              [](const vector<int>& part) { return part.size() <= 1; }),
                    partition.end());
    return partition;
}

vector<int> discriminate_code_set(const vector<vector<int>>& data) {
    int initial_col = choose_initial_column(data);
    vector<int> genes(data[0].size());
    for (int i = 0; i < genes.size(); ++i) {
        genes[i] = i;
    }

    vector<vector<int>> final_partition(2);
    for (int i = 0; i < data.size(); ++i) {
        final_partition[data[i][initial_col]].push_back(i);
    }

    vector<int> selected_cols = {initial_col};
    set<int> unselected_cols(genes.begin(), genes.end());
    unselected_cols.erase(initial_col);

    while (at_least_one(final_partition)) {
        int best_col = best_column(data, final_partition, unselected_cols);
        final_partition = split(final_partition, data, best_col);
        final_partition = remove_small_partition(final_partition);
        selected_cols.push_back(best_col);
        unselected_cols.erase(best_col);
    }

    int patient = final_partition[0][0];
    map<int, int> last_feature;
    for (int i = 0; i < data[0].size(); ++i) {
        if (data[patient][i] == 1) {
            int sum = 0;
            for (const auto& row : data) {
                sum += row[i];
            }
            last_feature[i] = sum;
        }
    }

    int best_feature = max_element(last_feature.begin(), last_feature.end(),
                                   [](const pair<int, int>& a, const pair<int, int>& b) {
                                       return a.second < b.second;
                                   })->first;
    selected_cols.push_back(best_feature);

    return selected_cols;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <data_file.csv>" << std::endl;
        return 1;
    }

    std::string data_file = argv[1];
    std::cout << "Processing file: " << data_file << std::endl;

    // Load data and gene list from the CSV file
    auto [data, gene_list] = load_csv(data_file);
    auto start_time = std::chrono::high_resolution_clock::now();
    // Assuming `discriminate_code_set` is defined elsewhere and returns the index set.
    std::vector<int> discos = discriminate_code_set(data);
    auto end_time = std::chrono::high_resolution_clock::now();

    std::string cancer_type = get_cancer_type(data_file);
//    std::string output_file = "DCS" + cancer_type + ".csv";

    std::vector<std::string> genes_result = get_genes_by_indices(gene_list, discos);
//    std::ofstream out_file(output_file);
//
//    if (!out_file.is_open()) {
//        std::cerr << "Error: Unable to open output file " << output_file << std::endl;
//        return 1;
//    }
//
//    // Write gene names to output file
//    for (const auto& gene : genes_result) {
//        out_file << gene << std::endl;
//    }
//
//    out_file.close();
        // Join the gene names with commas and output to the terminal
    std::cout << "Patients: " << data.size() << std::endl;
    std::cout << "Genes: " << data[0].size() << std::endl;
    std::cout << "DCS: ";
    for (size_t i = 0; i < genes_result.size(); ++i) {
        std::cout << genes_result[i];
        if (i != genes_result.size() - 1) {
            std::cout << ",";
        }
    }
    std::cout << std::endl;
    std::chrono::duration<double> elapsed_time = end_time - start_time;
    std::cout << "Time (seconds): " << elapsed_time.count() << std::endl;
    std::cout << "Size of DCS: " << discos.size() << std::endl;
    return 0;
}