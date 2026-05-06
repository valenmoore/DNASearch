#include <iostream>
#include <limits>
#include <chrono>
#include <map>
#include "json.hpp"
#include "FastaParser.h"
#include "Search.h"
#include "DNASequence.h"
#include "Species.h"


int main() {
    std::string fasta_path = "GCF_000005845.2";

    /*std::string query = "ACGATCTCGGCTCACTGCAAGGTCCGCCTCCCGGGTTCACGCCATTCTCCTGCCTCAGTCTCCCGAGTAGCTGGGACCACAGGCGCCCGCCACCATGCCCAGCTAGTTTTTTGTATTTTTGGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTGGATCTCCTGACCTCGTGATCCACCAGCCTCGGCCTCCCAAAGTGCTGGGATTACAGACGTGAGCCACCGTGCCCAGCTGAGAAAATGGGGTTTTCTAAATATACAATCATGTCATCTGCAAACAGAGACCATTTGACTTCCTCTCTTCCTATTTGAATACCCTTTATTTCTTTCTCTTGCCTCA";
    std::cout << std::endl;
    DNASequence querySeq = DNASequence(query);
    std::cout << "building search..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    Search s = Search(fasta_path);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "build duration: " << duration.count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    std::vector<int> indices = s.dumbSearch(s.getGenome(), querySeq, 0);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;

    std::cout << "dumb search duration: " << duration.count() << std::endl;
    for (int i : indices) {
        std::cout << i << std::endl;
    }
    std::cout << "-----------" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::unordered_set<int> things = s.smartSearch(querySeq, 0);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "smart search duration: " << duration.count() << std::endl;
    for (int i : things) {
        std::cout << i << std::endl;
    }*/

    std::vector<Species*> species;
    species.push_back(new Species("e coli", "GCF_000005845.2"));
    species.push_back(new Species("human", "NC_000001.11"));

    while (true) {
        std::cout << "\nEnter JSON query (or 'quit'): " << std::endl;
        std::cout << "Example: {\"query\": \"ATCG\", \"genomes\": [\"NC_000001.11\"], \"maxDist\": 2}" << std::endl;

        std::string line;
        std::getline(std::cin, line);
        if (line == "quit") break;

        // parse
        nlohmann::json body;
        try {
            body = nlohmann::json::parse(line);
        } catch (const nlohmann::json::parse_error& e) {
            std::cout << "Invalid JSON: " << e.what() << std::endl;
            continue;
        }

        // validate query
        if (!body.contains("query") || !body["query"].is_string()) {
            std::cout << "Missing or invalid 'query' field." << std::endl;
            continue;
        }
        std::string query = body["query"];

        // validate DNA
        bool valid = true;
        for (char c : query) {
            if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
                std::cout << "Invalid character '" << c << "' in sequence." << std::endl;
                valid = false;
                break;
            }
        }
        if (!valid) continue;

        // validate maxDist
        int maxDist = 3; // default
        if (body.contains("maxDist")) {
            if (!body["maxDist"].is_number_integer() || body["maxDist"] < 0) {
                std::cout << "Invalid 'maxDist', must be a non-negative integer." << std::endl;
                continue;
            }
            maxDist = body["maxDist"];
        }

        // validate genomes
        std::vector<Species*> toSearch;
        if (!body.contains("genomes") || body["genomes"] == "all") {
            for (const auto& s : species) toSearch.push_back(s);
        } else if (body["genomes"].is_array()) {
            for (const auto& g : body["genomes"]) {
                if (!g.is_string()) {
                    std::cout << "Genome ids must be strings." << std::endl;
                    continue;
                }
                std::string genomeId = g;
                for (const auto& s : species) {
                    if (s->getNcbiId() == genomeId) toSearch.push_back(s);
                }
            }
        } else {
            std::cout << "Invalid 'genomes' field, must be an array or 'all'." << std::endl;
            continue;
        }

        if (toSearch.empty()) {
            std::cout << "No valid genomes selected." << std::endl;
            continue;
        }

        // search and output as JSON
        nlohmann::json results;
        results["query"] = query;
        results["maxDist"] = maxDist;

        DNASequence querySeq(query);
        for (const auto& spec : toSearch) {
            auto start = std::chrono::high_resolution_clock::now();
            auto hits = spec->search(querySeq, maxDist);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> duration = end - start;

            results["results"][spec->getName()]["hits"] = std::vector<int>(hits.begin(), hits.end());
            results["results"][spec->getName()]["count"] = hits.size();
            results["results"][spec->getName()]["durationMs"] = duration.count();
        }

        std::cout << results.dump(2) << std::endl;
    }

    return 0;
}
