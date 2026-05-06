//
// Created by Valen Moore on 5/5/26.
//

#ifndef DNASEARCH_SPECIES_H
#define DNASEARCH_SPECIES_H
#include <string>
#include <utility>

#include "Search.h"

/**
 * Class to hold species name, id, and searcher for user interaction
 */
class Species {
    std::string name;
    std::string ncbiId;
    Search* searcher;

public:
    /**
     * Construct a species
     * @param name name of the species
     * @param ncbiId ncbiId of the species
     */
    Species(std::string name, const std::string &ncbiId) : name(std::move(name)), ncbiId(ncbiId) {
        searcher = new Search(ncbiId);
    }

    /**
     * Gets the name
     * @return species name
     */
    const std::string &getName() const { return name; }

    /**
     * Gets the NCBI id
     * @return id
     */
    const std::string &getNcbiId() const { return ncbiId; }

    /**
     * Searches for a query (wrapper function for search using species' genome
     * @param query the query to search for
     * @param maxDist maximum allowed hamming distance between match and query
     * @return a set of indices of matches found in the genome
     */
    [[nodiscard]] std::unordered_set<int> search(const DNASequence& query, int maxDist) const {
        return searcher->smartSearch(query, maxDist);
    }

    ~Species() {
        delete searcher;
    }
};


#endif //DNASEARCH_SPECIES_H