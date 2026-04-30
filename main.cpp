#include <iostream>
#include <limits>
#include <chrono>
#include "FastaParser.h"
#include "Search.h"
#include "DNASequence.h"

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
    std::string seq = FastaParser::readSequenceAsStr(fasta_path, 0, 248956422);
    // std::cout << seq << std::endl;
    std::string query = "TAGGAGGCAGAGCTGTCT";
    std::cout << "sequence loaded." << std::endl;
    DNASequence smartSeq = DNASequence(seq);
    std::cout << std::endl;
    DNASequence querySeq = DNASequence(query);
    std::cout << "building search...";
    auto start = std::chrono::high_resolution_clock::now();
    Search s = Search(smartSeq);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "build duration: " << duration.count() << std::endl;
    // smartSeq.print();
    querySeq.print();
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> indices = s.dumbSearch(smartSeq, querySeq, 2);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;

    std::cout << "dumb search duration: " << duration.count() << std::endl;
    for (int i : indices) {
        std::cout << i << std::endl;
    }
    std::cout << "-----------" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::unordered_set<int> things = s.smartSearch(querySeq, 2);
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "dumb search duration: " << duration.count() << std::endl;
    for (int i : things) {
        std::cout << i << std::endl;
    }
    /*Search s(seq);
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
    }*/

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
