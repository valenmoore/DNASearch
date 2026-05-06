//
// Created by Valen Moore on 5/5/26.
//

#ifndef DNASEARCH_SPECIES_H
#define DNASEARCH_SPECIES_H
#include <string>
#include <utility>

#include "Search.h"


class Species {
    std::string name;
    std::string ncbiId;
    Search* searcher;

public:
    Species(std::string name, const std::string &ncbiId) : name(std::move(name)), ncbiId(ncbiId) {
        searcher = new Search(ncbiId);
    }

    const std::string &getName() const { return name; }
    const std::string &getNcbiId() const { return ncbiId; }

    [[nodiscard]] std::unordered_set<int> search(const DNASequence& query, int maxDist) const {
        return searcher->smartSearch(query, maxDist);
    }

    ~Species() {
        delete searcher;
    }
};


#endif //DNASEARCH_SPECIES_H