//
// Created by Valen Moore on 4/16/26.
//

#include "Search.h"

#include <iostream>

std::vector<int> Search::dumbSearch(const std::string &query, int maxDist) {
    std::vector<int> out;
    for (int i = 0; i <= genomeSequence.length() - query.length(); i++) {
        std::string substr = genomeSequence.substr(i, query.length());
        if (hammingDistance(substr, query, maxDist) <= maxDist) {
            out.push_back(i);
        }
    }
    return out;
}

std::unordered_set<int> Search::smartSearch(const std::string &query, int maxDist) {
    std::vector<uint64_t> querySignature = computeSignature(query);
    std::unordered_set<int> output;
    for (const Tile& t : genomeIndex) {
        if (estimateSimilarity(t.signature, querySignature) > 0.1) {
            std::string tileSeq = genomeSequence.substr(t.position, 100);

            for (int i = 0; i <= tileSeq.length() - query.length(); i++) {
                if (hammingDistance(tileSeq.substr(i, query.length()), query, maxDist) <= maxDist) {
                    output.insert(t.position + i);
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

uint64_t Search::hashKmer(const std::string& kmer, const std::string& seed) {
    return std::hash<std::string_view>{}(kmer) ^ std::hash<std::string>{}(seed);
}

double Search::estimateSimilarity(const std::vector<uint64_t>& sig1,const std::vector<uint64_t>& sig2) {
    int matches = 0;
    for (int i = 0; i < numHashes; i++) {
        if (sig1[i] == sig2[i]) {
            matches++;
        }
    }
    return static_cast<double>(matches) / numHashes;
}

std::vector<uint64_t> Search::computeSignature(const std::string &sequence) {
    std::vector<uint64_t> signature(numHashes);

    for (int i = 0; i < numHashes; i++) {
        uint64_t minHash = UINT64_MAX;
        for (int j = 0; j <= sequence.length() - kLength; j++) {
            std::string_view kmer(sequence.data() + j, kLength);
            uint64_t hash = hashKmer(std::string(kmer), hashSeeds[i]);
            if (hash < minHash) minHash = hash;
        }
        signature[i] = minHash;
    }
    return signature;
}

int Search::hammingDistance(const std::string &subSeq, const std::string &query, int maxDist) {
    int dist = 0;
    for (int i = 0; i < subSeq.length(); i++) {
        if (subSeq[i] != query[i]) dist++;
        if (dist > maxDist) return dist; // prevent having to iterate further than necessary
    }
    return dist;
}
