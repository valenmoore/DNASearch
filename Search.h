//
// Created by Valen Moore on 4/16/26.
//

#ifndef DNASEARCH_SEARCH_H
#define DNASEARCH_SEARCH_H
#include <filesystem>
#include <iostream>
#include <random>
#include <utility>
#include <vector>
#include <string>
#include <unordered_set>

#include "DNASequence.h"


struct Tile {
    int position;
    std::vector<uint64_t> signature;
};

class Search {
private:
    const std::string genomeId;
    std::vector<uint64_t> hashSeeds;
    std::vector<Tile> genomeIndex;
    DNASequence genome;
    int numHashes = 100;
    int kLength = 11;
    int windowSize = 400;
    int stepSize = 200;

    std::vector<std::uint64_t> generateHashSeeds(int count) {
        std::vector<std::uint64_t> seeds;
        std::mt19937_64 gen(42);
        std::uniform_int_distribution<uint64_t> dis;

        for (int i = 0; i < count; i++) {
            seeds.push_back(dis(gen));
        }
        return seeds;
    }

    uint64_t hashKmer(uint64_t x, uint64_t seed);
    std::unordered_set<std::string> getKmers(const std::string& sequence);
    double estimateSimilarity(const std::vector<uint64_t>& sig1, const std::vector<uint64_t>& sig2);

    void buildIndex();
    void loadIndex(const std::string& path);
    void saveIndex(const std::string& savePath);
public:
    Search(DNASequence genome) : genome(std::move(genome)) {
        hashSeeds = generateHashSeeds(numHashes);
        buildIndex();
    }

    Search(const std::string& genomeId) : genome(FastaParser::readSequenceAsStr(genomeId)), genomeId(genomeId) {
        std::string filePath = "../cachedGenomes/" + genomeId + ".bin";
        if (std::filesystem::exists(filePath)) {
            std::cout << "File exists: " << filePath << std::endl;
            loadIndex(filePath);
        } else {
            std::cout << "could not find cache, building index" << std::endl;
            hashSeeds = generateHashSeeds(numHashes);
            buildIndex();
            saveIndex(filePath);
        }
    }
    std::vector<uint64_t> computeSignature(const DNASequence& sequence);
    std::vector<int> dumbSearch(const DNASequence& seq, const DNASequence& query, int maxDist=3);
    std::unordered_set<int> smartSearch(const DNASequence& query, int maxDist=3);

    std::vector<uint32_t> search(const std::string& sequence, const std::string& query, int maxDist=3);

    const DNASequence& getGenome() const { return genome; }
};


#endif //DNASEARCH_SEARCH_H