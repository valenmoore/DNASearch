//
// Created by Valen Moore on 4/20/26.
//

#ifndef DNASEARCH_DNASEQUENCE_H
#define DNASEARCH_DNASEQUENCE_H

#include <vector>

#include "FastaParser.h"


class DNASequence {
private:
    const std::vector<uint64_t> sequence;
    size_t numBases;
public:
    DNASequence(std::string dna) : sequence(FastaParser::dnaToBits(dna)) {
        numBases = dna.length();
    }
    DNASequence(uint64_t dna, size_t numBases) : sequence({dna}), numBases(numBases) {}
    DNASequence(std::vector<uint64_t> dna, size_t numBases) : sequence(std::move(dna)), numBases(numBases) {}

    [[nodiscard]] uint8_t getBase(int index) const;
    [[nodiscard]] uint64_t getSubseq(int start, int length) const;

    [[nodiscard]] size_t getNumBases() const { return numBases; };
    void print() const;
};


#endif //DNASEARCH_DNASEQUENCE_H