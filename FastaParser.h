//
// Created by Valen Moore on 4/15/26.
//

#ifndef DNASEARCH_FASTAPARSER_H
#define DNASEARCH_FASTAPARSER_H
#include <string>
#include <vector>


class FastaParser {
public:
    /**
     * Reads an NCBI genome sequence as bits, returning a vector of integers containing bits for the sequence
     * @param id the NCBI id of the genome
     * @return the vector of integers
     */
    static std::vector<uint64_t> readSequenceAsBits(const std::string& id);

    /**
     * Reads an NCBI genome sequence as a string of bases
     * @param id the NCBI id of the genome
     * @return the string of bases
     */
    static std::string readSequenceAsStr(const std::string& id);

    /**
     * Converts a string of DNA into a vector of bits
     * @param dna the string of bases
     * @return an array of ints containing the bits representing the DNA
     */
    static std::vector<uint64_t> dnaToBits(std::string& dna);

    /**
     * Converts a 32-base or less sequence of bits into a string
     * @param bits the sequence, in bits
     * @return the string version
     */
    static std::string bitsToDna(uint64_t bits);
private:
    static std::string &cleanSequence(std::string &sequence);
};


#endif //DNASEARCH_FASTAPARSER_H