//
// Created by Valen Moore on 4/16/26.
//

#include "Search.h"

#include <iostream>

std::vector<int> Search::dumbSearch(const DNASequence& seq, const DNASequence& query, int maxDist) {
    std::vector<int> out;
    int distance = 0;
    for (int i = 0; i <= seq.getNumBases() - query.getNumBases(); i++) {
        for (int j = 0; j < query.getNumBases(); j++) {
            if (seq.getBase(i + j) != query.getBase(j)) {
                distance++;
            }
            if (distance > maxDist) break;
        }
        if (distance <= maxDist) out.push_back(i);
        distance = 0;
    }
    return out;
}

void Search::buildIndex() {
    uint64_t mask = (1ULL << (kLength * 2)) - 1; // mask with size of kmer
    std::vector<uint64_t> tileSignature(numHashes);
    for (int i = 0; i <= genome.getNumBases() - windowSize; i += stepSize) {
        std::fill(tileSignature.begin(), tileSignature.end(), UINT64_MAX);
        uint64_t kmer = genome.getSubseq(i, kLength);
        for (int j = 0; j < windowSize; j++) {
            int baseIndex = j + i;
            kmer = (kmer << 2 | genome.getBase(baseIndex)) & mask; // shift over and insert new base to prevent recalculating kmer every time
            if (j >= kLength - 1) {
                uint64_t primaryHash = hashKmer(kmer, 0);

                for (int h = 0; h < numHashes; h++) {
                    uint64_t derivedHash = primaryHash ^ hashSeeds[h];
                    if (derivedHash < tileSignature[h]) {
                        tileSignature[h] = derivedHash;
                    }
                }
            }
        }
        Tile t;
        t.position = i;
        t.signature = tileSignature;
        genomeIndex.push_back(t);
    }
}

std::unordered_set<int> Search::smartSearch(const DNASequence& query, int maxDist) {
    std::vector<uint64_t> querySignature = computeSignature(query);
    std::unordered_set<int> output;
    for (const Tile& t : genomeIndex) {
        if (estimateSimilarity(t.signature, querySignature) > 0.01) {
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
                    output.insert(i);
                }
            }
        }
    }
    return output;
}

std::unordered_set<std::string> Search::getKmers(const std::string &sequence) {
    std::unordered_set<std::string> kmers;
    for (int i = 0; i <= sequence.length() - kLength; i++) {
        kmers.insert(sequence.substr(i, kLength));
    }
    return kmers;
}

uint64_t Search::hashKmer(uint64_t x, const uint64_t seed) {
    x ^= seed;

    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdLLU;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53LLU;
    x ^= x >> 33;

    return x;
}

double Search::estimateSimilarity(const std::vector<uint64_t>& sig1, const std::vector<uint64_t>& sig2) {
    int matches = 0;
    for (int i = 0; i < numHashes; i++) {
        if (sig1[i] == sig2[i]) {
            matches++;
        }
    }
    return static_cast<double>(matches) / numHashes;
}

std::vector<uint64_t> Search::computeSignature(const DNASequence& sequence) {
    std::vector<uint64_t> signature(numHashes, UINT64_MAX);
    uint64_t mask = (1ULL << (kLength * 2)) - 1;
    uint64_t kmer = 0;

    for (int i = 0; i < sequence.getNumBases(); i++) {
        kmer = ((kmer << 2) | sequence.getBase(i)) & mask;

        if (i >= kLength - 1) {
            uint64_t primaryHash = hashKmer(kmer, 0);

            for (int h = 0; h < numHashes; h++) {
                uint64_t derivedHash = primaryHash ^ hashSeeds[h];
                if (derivedHash < signature[h]) {
                    signature[h] = derivedHash;
                }
            }
        }
    }
    return signature;
}

