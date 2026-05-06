# Genome Search API

## Overview

This project using Locality Sensitive Hashing to query genomes for similar sequences of DNA. It is sort of reminiscent of NCBI's BLAST tool, but implemented from scratch in C++ (with many simplifications).

Rather than linearly searching the entire genome for exact matches, a genome is converted into an "index". The genome is tiled, and each tile is further split into kmers. A tile is represented by an array of hashes called the signature. Each hash is the minimum of 100 different permutations. 

For querying, the signature of the query is calculated, and then then each tile in the index is compared to the query to find the Jaccard Similarity between the two. If the similarity is above a threshold, the tile is searched for sequences that match the query within a given tolerance.

This technique is called [MinHash](https://en.wikipedia.org/wiki/MinHash), or Min-wise Independent Permutations Locality Sensitive Hashing, and it was proposed by Andrei Broder. It allows for simple integer comparisons to filter out dissimilar tiles, increasing efficiency on large datasets like genomes. I implemented a simplified version of MinHash for this project.

In addition to MinHash, I increased memory efficiency by storing DNA as bits instead of strings. Each base (A, T, C, or G) was represented by a two-bit sequence (00, 01, 10, or 11). Because one character in a string takes 8 bits, this method decreases memory usage by 4x by using only 2 bits / base. Operations with binary, like hashing, were also more time efficient than string manipulation.

## Usage
To use the project, run the main.cpp file. It will prompt you for a query in JSON format, which you can enter. The format is as follows:
```{"query": "ATCG", "genomes": ["NC_000001.11"], "maxDist": 2}```. The query represents the DNA query. The genomes represents the genomes to search. Replace the array with "all" to search all genomes in the cachedGenomes folder. The maxDist is the maximum allowed distance between the query and a potential match. Please note that the LSH parameters are optimized for long searches, of 50+ bases. Short queries would require architectural changes, but the main goal of this project was to search for and identify long strands of DNA.

To query new genomes, first locate the genome's NCBI ID. You can find the genome search [here](https://www.ncbi.nlm.nih.gov/home/genomes/); the ID is next to the "NCBI RefSeq assembly" label on the genome's main page. For example, the [human coronavirus genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000853505.1/) has an ID of GCF_000853505.1. Add this ID the array of species in main.cpp, along with the name of the organism. Running the project will automatically download and build the index for this new species. Building indices is slow, so it may take up to a minute for larger genomes.

Each query will return a JSON-formatted output (or an error message). The format will be as follows:

```
{
  "maxDist": 2,
  "query": "ACTCGG..."
  "results": {
    "E. Coli": {
      "count": 0,
      "durationMs": 1084.448625,
      "hits": []
    },
    "Homo Sapiens Chromosome 21": {
      "count": 1,
      "durationMs": 2327.06925,
      "hits": [
        230456712
      ]
    }
  }
}
```

where the results for each species/genome will be stored under that species/genome's name. The count is the number of matches; the duration is the time duration of the query, in milliseconds, and the array of "hits" is the locations of matches found on the chromosome.
