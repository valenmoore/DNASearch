//
// Created by Valen Moore on 4/15/26.
//

#ifndef DNASEARCH_FASTAPARSER_H
#define DNASEARCH_FASTAPARSER_H
#include <string>
#include <vector>


class FastaParser {
public:
    static std::vector<uint64_t> readSequenceAsBits(const std::string& id);
    static std::string readSequenceAsStr(const std::string& id);

    static std::vector<uint64_t> dnaToBits(std::string& dna);
    static std::string bitsToDna(uint64_t bits);
private:
    static std::string &cleanSequence(std::string &sequence);
};


#endif //DNASEARCH_FASTAPARSER_H