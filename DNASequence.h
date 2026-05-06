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
    /**
     * Constructs a DNA sequence object
     * @param dna the string of DNA
     */
    DNASequence(std::string dna) : sequence(FastaParser::dnaToBits(dna)) {
        numBases = dna.length();
    }

    /**
     * Constructs a DNA sequence object
     * @param dna a vector of binary that represents the DNA
     * @param numBases the number of bases
     */
    DNASequence(std::vector<uint64_t> dna, size_t numBases) : sequence(std::move(dna)), numBases(numBases) {}

    /**
     * Gets the base at an index
     * @param index index of the base (same system as string indexing)
     * @return the base, as the first two bits in a uint8_t
     */
    [[nodiscard]] uint8_t getBase(int index) const;

    /**
     * Gets a subsequence of the DNA (cannot get subsequences longer than 32 bases)
     * @param start the start index
     * @param length length of the window
     * @return the subsequence, packed into one uint64_t as bits
     */
    [[nodiscard]] uint64_t getSubseq(int start, int length) const;

    /**
     * Gets the number of bases in the sequence
     * @return the number of bases
     */
    [[nodiscard]] size_t getNumBases() const { return numBases; };

    /**
     * Prints the sequence
     */
    void print() const;
};


#endif //DNASEARCH_DNASEQUENCE_H