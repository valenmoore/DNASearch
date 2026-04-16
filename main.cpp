#include <iostream>
#include <bitset>
#include <chrono>
#include "FastaParser.h"
#include "Search.h"

double estimateSimilarity(const std::vector<uint64_t>& sig1,
                         const std::vector<uint64_t>& sig2) {
    int matches = 0;
    for (int i = 0; i < 100; i++) {
        if (sig1[i] == sig2[i]) {
            matches++;
        }
    }
    return static_cast<double>(matches) / 100;
}

int main() {
    std::string fasta_path = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
    std::string seq = FastaParser::readSequenceAsStr(fasta_path, 50000000, 10000000);
    std::string query = "TTCCCAAATGGGGCAGAAGA";
    std::cout << "sequence loaded." << std::endl;
    std::cout << "building search...";
    auto start = std::chrono::high_resolution_clock::now();
    Search s(seq);
    std::cout << "end search build" << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "build duration: " << duration.count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    auto results = s.smartSearch(query, 5);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "smart latency: " << duration.count() << std::endl;

    for (int r : results) {
        std::cout << r << std::endl;
        std::cout << seq.substr(r, query.length()) << std::endl;
    }


    std::cout << "--other other other dumb dumb way--" << std::endl;
    start = std::chrono::high_resolution_clock::now();

    auto otherResults = s.dumbSearch(query, 5);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "dumb latency: " << duration.count() << std::endl;

    for (int r : results) {
        std::cout << r << std::endl;
        std::cout << seq.substr(r, query.length()) << std::endl;
    }

    /*auto start = std::chrono::high_resolution_clock::now();
    std::vector<int> v = Search::dumbSearch(seq, query, 5);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    for (int i : v) {
        std::cout << "similar sequence found at " <<  i << std::endl;
        std::cout << "sequence: " << seq.substr(i, query.length()) << std::endl;
    }
    std::cout << "duration (millis) " << duration << std::endl;*/

    return 0;
}
