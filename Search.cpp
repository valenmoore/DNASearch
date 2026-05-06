//
// Created by Valen Moore on 4/16/26.
//

#include "Search.h"

#include <fstream>
#include <iostream>

std::vector<int> Search::dumbSearch(const DNASequence& seq, const DNASequence& query, int maxDist) {
    std::vector<int> out;
    int distance = 0;

    // compare each window of the sequence to the query until the distance surpasses the maximum
    for (int i = 0; i <= seq.getNumBases() - query.getNumBases(); i++) {
        for (int j = 0; j < query.getNumBases(); j++) {
            if (seq.getBase(i + j) != query.getBase(j)) {
                distance++;
            }
            if (distance > maxDist) break;
        }
        if (distance <= maxDist) out.push_back(i); // match found if distance is less than max
        distance = 0;
    }
    return out;
}

/**
 * Splits the genome into tiles and calculates the signature for each tile using min hashing
 */
void Search::buildIndex() {
    uint64_t mask = (1ULL << (kLength * 2)) - 1; // mask with size of kmer
    std::vector<uint64_t> tileSignature(numHashes);

    for (int i = 0; i <= genome.getNumBases() - windowSize; i += stepSize) {
        std::fill(tileSignature.begin(), tileSignature.end(), UINT64_MAX); // to fill with hash values
        uint64_t kmer = 0;
        for (int j = 0; j < windowSize; j++) {
            int baseIndex = j + i;
            if (baseIndex >= genome.getNumBases()) break;
            kmer = (kmer << 2 | genome.getBase(baseIndex)) & mask; // shift over and insert new base to prevent recalculating kmer every time
            if (j >= kLength - 1) {
                // if kmer is the right length, hash it numHashes times and then insert the minimum hash into the tileSignature
                uint64_t hash = 0;
                for (int h = 0; h < numHashes; h++) {
                    hash = hashKmer(kmer, hashSeeds[h]);
                    if (hash < tileSignature[h]) {
                        tileSignature[h] = hash;
                    }
                }
            }
        }

        // add to the genome index for use in smart search
        Tile t;
        t.position = i;
        t.signature = tileSignature;
        genomeIndex.push_back(t);
    }
    std::cout << "built" << std::endl;
}

std::unordered_set<int> Search::smartSearch(const DNASequence& query, int maxDist) {
    std::vector<uint64_t> querySignature = computeSignature(query);
    std::unordered_set<int> output;

    for (const Tile& t : genomeIndex) {
        // iterate over tiles and estimate similarity between tile and query
        if (estimateSimilarity(t.signature, querySignature) > 0.00) {
            // if tile is similar, execute a "dumb search" over the tile to find for exact matches
            int searchStart = t.position;
            int searchEnd = t.position + windowSize;

            if (searchEnd > (int)genome.getNumBases() - (int)query.getNumBases()) {
                searchEnd = (int)genome.getNumBases() - (int)query.getNumBases();
            }

            for (int i = searchStart; i < searchEnd; i++) {
                int distance = 0;
                for (int j = 0; j < (int)query.getNumBases(); j++) {
                    if (genome.getBase(i + j) != query.getBase(j)) {
                        distance++;
                    }
                    if (distance > maxDist) break;
                }
                if (distance <= maxDist) {
                    output.insert(i); // add to output if within hamming distance tolerance
                }
            }
        }
    }

    return output;
}

/**
 * Get all kmers (windows) of a sequence
 * @param sequence the sequence
 * @return the kmers
 */
std::unordered_set<std::string> Search::getKmers(const std::string &sequence) {
    std::unordered_set<std::string> kmers;
    for (int i = 0; i <= sequence.length() - kLength; i++) {
        kmers.insert(sequence.substr(i, kLength));
    }
    return kmers;
}

/**
 * Hashes a kmer (hash function from Claude)
 * @param x the kmer
 * @param seed the seed of the hash, from hashSeeds
 * @return a hash value as an integer
 */
uint64_t Search::hashKmer(uint64_t x, const uint64_t seed) {
    x ^= seed;

    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdLLU;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53LLU;
    x ^= x >> 33;

    return x;
}

/**
 * Estimates similarity between two signatures using Jaccard Similarity
 * @param sig1 the first signature of hashes
 * @param sig2 the second signature of hashes
 * @return the similarity as a decimal
 */
double Search::estimateSimilarity(const std::vector<uint64_t>& sig1, const std::vector<uint64_t>& sig2) {
    int matches = 0;
    for (int i = 0; i < numHashes; i++) {
        if (sig1[i] == sig2[i]) {
            matches++;
        }
    }
    return static_cast<double>(matches) / numHashes;
}

/**
 * Hashes a sequence to compute its signature
 * @param sequence the sequence
 * @return the signature, as an array of hashes
 */
std::vector<uint64_t> Search::computeSignature(const DNASequence& sequence) {
    std::vector<uint64_t> signature(numHashes, UINT64_MAX);
    uint64_t mask = (1ULL << (kLength * 2)) - 1;
    uint64_t kmer = 0;

    for (int i = 0; i < sequence.getNumBases(); i++) {
        kmer = ((kmer << 2) | sequence.getBase(i)) & mask;

        if (i >= kLength - 1) {
            uint64_t hash = 0;
            for (int h = 0; h < numHashes; h++) {
                hash = hashKmer(kmer, hashSeeds[h]);
                if (hash < signature[h]) {
                    signature[h] = hash;
                }
            }
        }
    }
    return signature;
}

void Search::loadIndex(const std::string& path) {
    std::ifstream inFile(path, std::ios::binary);

    int savedNumHashes;
    inFile.read(reinterpret_cast<char*>(&savedNumHashes), sizeof(savedNumHashes));

    numHashes = savedNumHashes;

    hashSeeds.resize(numHashes);
    for (int i = 0; i < numHashes; i++) {
        inFile.read(reinterpret_cast<char*>(&hashSeeds[i]), sizeof(hashSeeds[i]));
    }
    int tileCount = 0;
    while (inFile) {
        Tile t;
        if (!inFile.read(reinterpret_cast<char*>(&t.position), sizeof(t.position))) break;

        t.signature.resize(numHashes);
        inFile.read(reinterpret_cast<char*>(t.signature.data()), numHashes * sizeof(uint64_t));

        genomeIndex.push_back(t);
        tileCount++;
    }
}

void Search::saveIndex(const std::string& savePath) {
    std::cout << "saving at " << savePath << std::endl;
    std::ofstream outFile(savePath, std::ios::binary);
    if (outFile.is_open()) {
        outFile.write(reinterpret_cast<const char*>(&numHashes), sizeof(numHashes));

        for (const auto& seed : hashSeeds) {
            outFile.write(reinterpret_cast<const char*>(&seed), sizeof(seed));
        }

        for (const auto& t : genomeIndex) {
            outFile.write(reinterpret_cast<const char*>(&t.position), sizeof(t.position));
            outFile.write(reinterpret_cast<const char*>(t.signature.data()),
                          t.signature.size() * sizeof(uint64_t));
        }
        outFile.close();
    }
}
