#include <iostream>
#include <limits>
#include <chrono>
#include "FastaParser.h"
#include "Search.h"
#include "DNASequence.h"


int main() {
    std::string fasta_path = "NC_000001.11";

    std::string query = "ACGATCTCGGCTCACTGCAAGGTCCGCCTCCCGGGTTCACGCCATTCTCCTGCCTCAGTCTCCCGAGTAGCTGGGACCACAGGCGCCCGCCACCATGCCCAGCTAGTTTTTTGTATTTTTGGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTGGATCTCCTGACCTCGTGATCCACCAGCCTCGGCCTCCCAAAGTGCTGGGATTACAGACGTGAGCCACCGTGCCCAGCTGAGAAAATGGGGTTTTCTAAATATACAATCATGTCATCTGCAAACAGAGACCATTTGACTTCCTCTCTTCCTATTTGAATACCCTTTATTTCTTTCTCTTGCCTCA";
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
    }

    return 0;
}
