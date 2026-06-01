#include <iostream>
#include <string>
#include <numeric>
#include "toml.hpp" 
#include <vector>
#include <cmath>

// Strictly typed struct to hold the application state
struct BioLSHasherConfig {
    uint32_t Mutation_Model;
    uint32_t Mutation_Expression_Type;

    uint32_t Short_Sequence_Length;
    uint32_t Long_Sequence_Length;
    
    uint32_t Coll_Num_Hashes_To_And_Sweep_Start;
    uint32_t Coll_Num_Hashes_To_And_Sweep_End;
    uint32_t Coll_Num_Hashes_To_Or_Sweep_Start;
    uint32_t Coll_Num_Hashes_To_Or_Sweep_End;

    uint32_t Sequence_Length_For_Similarity_Search;
    uint32_t Total_Sequences_In_Database_To_Search_From;
    uint32_t Total_Query_Sequences_To_Search;
    uint32_t Number_of_Runs_to_get_Average;

    uint32_t SimSearch_Num_Hashes_To_And_Sweep_Start;
    uint32_t SimSearch_Num_Hashes_To_And_Sweep_End;
    uint32_t SimSearch_Num_Hashes_To_Or_Sweep_Start;
    uint32_t SimSearch_Num_Hashes_To_Or_Sweep_End;

    float Similarity_Threshold;

    bool are_Bases_Drawn_From_Uniform_Distribution;
    std::vector<double> Categorical_Distribution_Probabilities;
};



// Isolated parsing function
BioLSHasherConfig load_config(const std::string& filepath) {
    try {
        // Parse the file into a TOML table
        toml::table config = toml::parse_file(filepath);

        std::vector<double> extracted_categorical_probs;
        
        if (toml::array* arr = config["Common_Configuration"]["Categorical_Distribution_Probabilities"].as_array()) {
            extracted_categorical_probs.reserve(arr->size());
            for (toml::node& node : *arr) {
                // value_or(0.0) forces the compiler to deduce double, avoiding -Wnarrowing
                extracted_categorical_probs.push_back(node.value_or(0.25)); 
            }
            
            // checking if sum = 1.0
            double sum = std::accumulate(extracted_categorical_probs.begin(), extracted_categorical_probs.end(), 0.0);
            if (std::abs(sum - 1.0) > 1e-6) {
                throw std::domain_error("Probability distribution do not sum to 1.0");
            }
        } else {
            throw std::runtime_error("Configuration error: 'Probabilities' array is missing or malformed.");
        }

        // Extract values using type-safe accessors and fallback defaults
        return BioLSHasherConfig{
            config["Mutation_Model_Configuration"]["Mutation_Model"].value_or(1u),
            config["Mutation_Model_Configuration"]["Mutation_Expression_Type"].value_or(0u),
            config["Collision_Test_Configuration"]["Short_Sequence_Length"].value_or(50u),
            config["Collision_Test_Configuration"]["Long_Sequence_Length"].value_or(100u),
            config["Collision_Test_Configuration"]["Coll_Num_Hashes_To_And_Sweep_Start"].value_or(1u),
            config["Collision_Test_Configuration"]["Coll_Num_Hashes_To_And_Sweep_End"].value_or(1u),
            config["Collision_Test_Configuration"]["Coll_Num_Hashes_To_Or_Sweep_Start"].value_or(1u),
            config["Collision_Test_Configuration"]["Coll_Num_Hashes_To_Or_Sweep_End"].value_or(1u),
            config["Similarity_Search_Test_Configuration"]["Sequence_Length_For_Similarity_Search"].value_or(50u),
            config["Similarity_Search_Test_Configuration"]["Total_Sequences_In_Database_To_Search_From"].value_or(10000u),
            config["Similarity_Search_Test_Configuration"]["Total_Query_Sequences_To_Search"].value_or(100u),
            config["Similarity_Search_Test_Configuration"]["Number_of_Runs_to_get_Average"].value_or(5u),
            config["Similarity_Search_Test_Configuration"]["SimSearch_Num_Hashes_To_And_Sweep_Start"].value_or(1u),
            config["Similarity_Search_Test_Configuration"]["SimSearch_Num_Hashes_To_And_Sweep_End"].value_or(1u),
            config["Similarity_Search_Test_Configuration"]["SimSearch_Num_Hashes_To_Or_Sweep_Start"].value_or(1u),
            config["Similarity_Search_Test_Configuration"]["SimSearch_Num_Hashes_To_Or_Sweep_End"].value_or(1u),
            config["Similarity_Search_Test_Configuration"]["Similarity_Threshold"].value_or(0.9f),
            config["Common_Configuration"]["are_Bases_Drawn_From_Uniform_Distribution"].value_or(true),
            extracted_categorical_probs
        };
    } catch (const toml::parse_error& err) {
        std::cerr << "Configuration parsing failed at " 
                  << err.source().begin << ":\n" << err.description() << "\n";
        throw; 
    }
}

int main() {
    const std::string config_path = "config.toml";
    
    try {
        // 1. Load the configuration
        BioLSHasherConfig config = load_config(config_path);
        
        // 3. Output loaded state
        std::cout << "--- Benchmarking Configuration Loaded ---\n";
        std::cout << "Mutation Model: " << config.Mutation_Model << "\n";
        std::cout << "Mutation Expression Type: " << config.Mutation_Expression_Type << "\n";
        std::cout << "Short Sequence Length: " << config.Short_Sequence_Length << "\n";
        std::cout << "Long Sequence Length: " << config.Long_Sequence_Length << "\n";
        std::cout << "Collision AND Sweep Range: [" 
                  << config.Coll_Num_Hashes_To_And_Sweep_Start << ", " 
                  << config.Coll_Num_Hashes_To_And_Sweep_End << "]\n";
        std::cout << "Collision OR Sweep Range: [" 
                  << config.Coll_Num_Hashes_To_Or_Sweep_Start << ", " 
                  << config.Coll_Num_Hashes_To_Or_Sweep_End << "]\n";
        std::cout << "Similarity Search AND Sweep Range: [" 
                  << config.SimSearch_Num_Hashes_To_And_Sweep_Start << ", " 
                  << config.SimSearch_Num_Hashes_To_And_Sweep_End << "]\n";
        std::cout << "Similarity Search OR Sweep Range: [" 
                  << config.SimSearch_Num_Hashes_To_Or_Sweep_Start << ", " 
                  << config.SimSearch_Num_Hashes_To_Or_Sweep_End << "]\n";
        std::cout << "Similarity Threshold: " << config.Similarity_Threshold << "\n";
        std::cout << "Bases Drawn From Uniform Distribution: " 
                  << (config.are_Bases_Drawn_From_Uniform_Distribution ? "Yes" : "No") << "\n";
        // std::cout << "Categorical Distribution Probabilities: [";
        // for (size_t i = 0; i < config.categorical_distribution_probabilities.size(); ++i) {
        //     std::cout << config.categorical_distribution_probabilities[i];
        //     if (i < config.categorical_distribution_probabilities.size() - 1) {
        //         std::cout << ", ";
        //     }
        // }        std::cout << "]\n";
        // std::cout << "-----------------------------------------\n";
        
        // Core execution logic (e.g., collision curve analysis, c-ANN search) 
        // proceeds here using the populated 'config' struct.
        
    } catch (const std::exception& e) {
        std::cerr << "Fatal error during initialization: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}