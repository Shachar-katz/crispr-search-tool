#include "functions.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"

// step one functions bellow:

void findKmersInFile(MultiFormatFileReader& fileReader, unordered_map<string,int>& globalKmerMap, int seedK, int minK, int legitimateSpacer){
    // line variable temporerally holds the reads
    string line;
    int progressCounter = 0;
    // we are looping over every read
    while (fileReader.getNextLine(line)) {
        progressCounter++;
        if (progressCounter % 100000 == 0){
            cout << "Procession line: " << progressCounter << endl;
        }
        if (skipThisLine(line)){
            cerr << "faulty line" << endl;
            continue;
        }
        // every line we create an empty Smap that maps from an Smer to vect of indecies in the line.
        unordered_map<string,vector<int>> singleLineMapSeedKToIdx;
        // we also create an empty set of unique pre "vetted" Kmers in this line.
        set<string> uniqueKmersInLine;
        // populating Smap for this line.
        findSeedPattern(line, singleLineMapSeedKToIdx, seedK);
        // for every Smer in the map, if it appears more then once in the line:
        // we check if it can be expanded to a Kmer and if so that Kmer is added to the set of unique Kmers
        for (const auto& [Smer, idxs] : singleLineMapSeedKToIdx){
            if (singleLineMapSeedKToIdx[Smer].size() > 1){
                expandSeedToKmer(line, Smer, idxs, minK, uniqueKmersInLine, legitimateSpacer);
            }
        }
        // for every unique Kmer found on the line:
        // 1) we Determine the canonical form of the K-mer to ensure consistency in representation 
        //     * choosing between reverse complement and the Kmer to not add both to the global Kmer map
        // 2) then we increment the count of the one that was selected in the global map (or adding it for the first time) 
        for (const auto& uniqueKmer : uniqueKmersInLine) {
            string uniqueKmerOrReverse = pickKey(uniqueKmer);
            globalKmerMap[uniqueKmerOrReverse]++;
        }
    }
}

bool skipThisLine(const string& read){
    int A = 0;
    int T = 0;
    int C = 0;
    int G = 0;
    for (char c : read){
        switch(c){
            case 'A': A++; break;
            case 'T': T++; break;
            case 'C': C++; break;
            case 'G': G++; break;
            default: break;
        }
    }
    int countOfMostFrequentNuc = max({A, T, C, G});
    double ratio = (double)countOfMostFrequentNuc / (double)read.size();
    if (ratio >= 0.5) {
        // Skip expansions for this line
        return true; 
    }
    else{
        return false;
    }
}

// this function populates the Smap
void findSeedPattern(string line, unordered_map<string,vector<int>>& singleLineMapSeedKToIdx, int seedK){
    for (int j = 0; j < (line.size() - seedK); j++){
        string key = line.substr(j,seedK);
        if (key.size() < seedK){
            break;
        }
        singleLineMapSeedKToIdx[key].emplace_back(j);
    }
}

// this function recives start and end location of 2 Kmers in expansion and 
// varifies that they are not overlapping
bool notOverlapping(int idxStartPotential, int idxStartCompare, int idxEndPotential, int idxEndCompare, int legitimateSpacer)
{
    return ((idxEndPotential + legitimateSpacer) < idxStartCompare) || ((idxEndCompare + legitimateSpacer) < idxStartPotential);
}

// this function populates the unique Kmer set
void expandSeedToKmer(const string& line, string Smer, vector<int> SmerIdxVect , int minK, set<string>& uniqueKmersInLine,int legitimateSpacer){
    // Track unique K-mers for this particular Smer on this particular line
    set<string> UniqueKmersFromSmer; 
    // we iterate over the index vector of the appearences of this Smer, each index gets a turn to be the potential Kmer.
    for (int potentialKmer = 0; potentialKmer < SmerIdxVect.size(); potentialKmer++){
        // we initilize an empty "best fitting Kmer for this index of choice"
        string bestK = "";
        // then we iterate over all the other indecies comparing the Smer appearance there to our "Potential Kmer"
        for (int comparisionKmer = 0; comparisionKmer < SmerIdxVect.size(); comparisionKmer++ ){
            // we verify that we are no trying to expand the same location
            if (comparisionKmer == potentialKmer){
                continue;
            }
            // we set the indecies of the start and end at our potential Kmer and our comparision
            int idxStartPotential = SmerIdxVect.at(potentialKmer);
            int idxEndPotential = SmerIdxVect.at(potentialKmer) + Smer.length() - 1;
            
            int idxStartCompare = SmerIdxVect.at(comparisionKmer);
            int idxEndCompare = SmerIdxVect.at(comparisionKmer) + Smer.length() - 1;
            
            // we set flags = true to represent are the 2 occurances equal at the start and end.
            bool equalAtStart = true;
            bool equalAtEnd = true;
            // while the 2 locations are equal at the start or end and the indecies are not overlapping:
            while((equalAtStart || equalAtEnd) && notOverlapping(idxStartPotential,idxStartCompare,idxEndPotential,idxEndCompare,legitimateSpacer)){
                // as long as the indecies would valid at the start if we decremented and they are still equal at the start:
                if (equalAtStart && idxStartPotential > 0 && idxStartCompare > 0){
                    // if the "next" nucleotid is equal on both occurances we decrement the position and continue expending
                    if (line[idxStartPotential - 1] == line[idxStartCompare - 1]){
                        idxStartPotential-- ;
                        idxStartCompare-- ;
                    }
                    // else we turn the equal at start flag to false.
                    else{
                        equalAtStart = false;
                    }
                }
                else{
                    equalAtStart = false;
                }
                // after one end expansion if they overlapp we want to stop
                if(!notOverlapping(idxStartPotential,idxStartCompare,idxEndPotential,idxEndCompare,legitimateSpacer)) {break;}
                // as long as the indecies would be valid at the end if we incremented and they are still equal at the end:
                if (equalAtEnd && (idxEndPotential + 1) < line.length() && (idxEndCompare + 1) < line.length()){
                    // if the "next" nucleotid is equal on both occurances we increment the position and continue expending
                    if (line[idxEndPotential + 1] == line[idxEndCompare + 1]){
                        idxEndPotential++ ;
                        idxEndCompare++ ;
                    }
                    // else we turn the equal at end flag to false.
                    else{
                        equalAtEnd = false;
                    }
                }
                else{
                    equalAtEnd = false;
                }
                
            } 
            // when we exist this loop we know that either:
            // the nucleotids dont match anymore or, they have reached a point that if expanded more they would overlap

            // once both flags are false and we are done expanding to both sides we can generate a Kmer:
            // we get the length of the (still "potential") Kmer
            int KmerLen = idxEndPotential - idxStartPotential + 1;

            // as long as the indecies are valid we generate the substring Kmer 
            if (idxStartPotential >= 0 && idxEndPotential < line.size()) {
                string Kmer = line.substr(idxStartPotential, KmerLen);
                // if this potential Kmer follows requirements :
                // (is indeed a Kmer and is longer then the best Kmer so far)
                // we make it the best Kmer
                if(KmerLen >= minK && KmerLen > bestK.length()){
                    bestK = Kmer;
                }
            }
        }
        // once we have compared every instance of the "potential Kmer" Smer:
        // if we got a Best Kmer we add it to the unique Kmers from Smer Set
        if(bestK.size() > 0){
            UniqueKmersFromSmer.insert(bestK);
        }
            
    }
    // once we have compared every Smer to all the other Smers:
    // We add the best match for each instance (a Kmer) to the unique Kmers in line Set.
    for (const auto& uniqueKmer : UniqueKmersFromSmer) {
        uniqueKmersInLine.insert(uniqueKmer);
    }
}
