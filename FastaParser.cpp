//
// Created by Valen Moore on 4/15/26.
//

#include "FastaParser.h"
#include <cpr/cpr.h>
#include <climits>
#include <iostream>

uint64_t FastaParser::readSequenceAsBits(std::string &fileName, int start, int windowLength) {
    std::string str = readSequenceAsStr(fileName, start, windowLength);
    return dnaToBits(str);
}

std::string FastaParser::readSequenceAsStr(std::string &fileName, int start, int windowLength) {
    cpr::Response r = cpr::Get(cpr::Url{fileName},
                           cpr::Parameters{{"db", "nuccore"},
                                           {"id", "NC_000022.11"},
                                           {"seq_start", std::to_string(start)},
                                           {"seq_stop", std::to_string(start + windowLength)},
                                           {"rettype", "fasta"}});

    if (r.status_code != 200) {
        throw std::invalid_argument("Error reading sequence");
    }
    cleanSequence(r.text);
    return r.text;
}

uint64_t FastaParser::dnaToBits(std::string& dna) {
    uint64_t bits = 0;

    for (char c : dna) {
        bits <<= 2; // shift over to make space
        switch (toupper(c)) {
            case 'A':
                bits |= 0; // binary 00
                break;
            case 'C':
                bits |= 1; // binary 01
                break;
            case 'G':
                bits |= 2; // binary 10
                break;
            case 'T':
                bits |= 3; // binary 11
                break;
            default:
                throw std::invalid_argument("Invalid character found in dnaToBits");
        }
    }
    return bits;
}

std::string FastaParser::bitsToDna(uint64_t bits) {
    return "";
}

std::string &FastaParser::cleanSequence(std::string &sequence) {
    sequence.erase(0, sequence.find('\n') + 1); // remove line with >sequence.erase(std::remove(sequence.begin(), sequence.end(), '\n'), sequence.end());
    sequence.erase(std::remove(sequence.begin(), sequence.end(), '\n'), sequence.end());
    sequence.erase(std::remove(sequence.begin(), sequence.end(), '\r'), sequence.end()); // handle windows line endings too
    sequence.erase(std::remove(sequence.begin(), sequence.end(), 'N'), sequence.end()); // Ns are telomeres which I don't care about
    return sequence;
}
