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

/**
 * A struct to hold a tile's position and the array of hashes that make up its signature
 */
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
    std::vector<uint64_t> computeSignature(const DNASequence& sequence);

    void buildIndex();
    void loadIndex(const std::string& path);
    void saveIndex(const std::string& savePath);
public:
    /**
     * Initialize a search class
     * @param genomeId the NCBI id of the genome, e.g. NC_000001.11
     */
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

    /**
     * A O(N*M)-ish search that iterates over every possible window of the genome; should rarely be used instead of smart search
     * @param seq input sequence to iterate over
     * @param query the query sequence
     * @param maxDist the maximum hamming distance to query for
     * @return the indices of each of the matches in the genome
     */
    std::vector<int> dumbSearch(const DNASequence& seq, const DNASequence& query, int maxDist=3);

    /**
     * Uses Locality Sensitive Hashing/clustering to search over a genome efficiently
     * @param query the query sequence
     * @param maxDist the maximum hamming distance to query for
     * @return the indices of each of the matches in the genome
     */
    std::unordered_set<int> smartSearch(const DNASequence& query, int maxDist=3);

    /**
     * Gets the genome sequence
     * @return the genome sequence as a DNASequence
     */
    const DNASequence& getGenome() const { return genome; }
};


#endif //DNASEARCH_SEARCH_H