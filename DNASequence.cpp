//
// Created by Valen Moore on 4/20/26.
//

#include "DNASequence.h"
#include <string>
#include <utility>
#include <vector>
#include <iostream>

uint8_t DNASequence::getBase(const int index) const {
    if (index >= numBases) throw std::runtime_error("Sequence index out of bounds");
    int bucketIndex = index / 32;
    int bitPos = 62 - ((index % 32) * 2);
    uint64_t bucket = sequence[bucketIndex];
    return bucket >> bitPos & 3; // shift everything over by the index and mask out two characters
}

uint64_t DNASequence::getSubseq(int start, int length) const {
    if (start + length >= numBases) throw std::runtime_error("Sequence index out of bounds");
    uint64_t result = 0;
    for (int i = start; i < start + length; i++) {
        result <<= 2; // shift over two spaces for new character
        result |= getBase(i);
    }
    return result;
}

void DNASequence::print() const {
    for (int i = 0; i < numBases; i++) {
        std::cout << std::bitset<2>(getBase(i)) << " ";
    }
    std::cout << std::endl;
}
