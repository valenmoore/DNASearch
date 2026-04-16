//
// Created by Valen Moore on 4/16/26.
//

#ifndef DNASEARCH_SEARCH_H
#define DNASEARCH_SEARCH_H
#include <random>
#include <vector>
#include <string>
#include <unordered_set>


struct Tile {
    int position;
    std::vector<uint64_t> signature;
};

class Search {
private:
    std::vector<std::string> hashSeeds;
    std::vector<Tile> genomeIndex;
    std::string genomeSequence;
    int numHashes = 100;
    int kLength = 4;
    int windowSize = 100;
    int stepSize = 50;

    static int hammingDistance(const std::string& subSeq, const std::string& query, int maxDist);

    std::vector<std::string> generateHashSeeds(int count) {
        std::vector<std::string> seeds;
        std::mt19937 gen(42);
        std::uniform_int_distribution<> dis(0, 3);
        const char bases[] = {'A', 'T', 'C', 'G'};

        for (int i = 0; i < count; i++) {
            std::string seed = "";
            for (int j = 0; j < 20; j++) {
                seed += bases[dis(gen)];
            }
            seeds.push_back(seed);
        }
        return seeds;
    }

    uint64_t hashKmer(const std::string& kmer, const std::string& seed);
    std::unordered_set<std::string> getKmers(const std::string& sequence);
    double estimateSimilarity(const std::vector<uint64_t>& sig1, const std::vector<uint64_t>& sig2);
    void buildIndex(const std::string& chromosome) {
        for (int i = 0; i <= chromosome.length() - windowSize; i += stepSize) {
            std::string window = chromosome.substr(i, windowSize);
            genomeIndex.push_back({i, computeSignature(window)});
        }
    }
public:
    Search(const std::string& sequence) {
        genomeSequence = sequence;
        hashSeeds = generateHashSeeds(numHashes);
        buildIndex(genomeSequence);
    }
    std::vector<uint64_t> computeSignature(const std::string& sequence);
    std::vector<int> dumbSearch(const std::string& query, int maxDist=3);
    std::unordered_set<int> smartSearch(const std::string& query, int maxDist=3);

    std::vector<uint32_t> search(const std::string& sequence, const std::string& query, int maxDist=3);
};


#endif //DNASEARCH_SEARCH_H