#include "functions.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"

// step one functions bellow:

void findKmersInFile(MultiFormatFileReader& fileReader, 
                     unordered_map<string,int>& globalKmerMap, 
                     int seedK, 
                     int minK, 
                     int legitimateSpacer,
                     unordered_map<string,double>& stats, 
                     bool strict, 
                     bool preStrict, 
                     ofstream& logFile,
                     int interval,
                     int maxK)
    {
    // line variable temporerally holds the reads
    string line;
    // statistics Vars
    int progressCounter = 0;
    int numReadsWithRepeats = 0;
    int faultyLine = 0;
    // we are looping over every read
    while (fileReader.getNextLine(line)) {
        // statistics and progress managment:
        progressCounter++;
        if (progressCounter % interval == 0){
            cout << "Procession line: " << progressCounter << endl;
            logFile << "Procession line: " << progressCounter << endl;
            // also output faulty lines
            logFile << faultyLine << " faulty lines" << endl;
            cerr << faultyLine << " faulty lines" << endl;
            faultyLine = 0;
        }
        // check for faulty lines
        if ((preStrict == 1 && skipThisLine(line, 0.5)) || line.length() <= (2 * minK + legitimateSpacer + 2)){
            faultyLine++;
            continue;
        }
        
        // every line we create an empty smap that maps from an smer to vect of indecies in the line.
        unordered_map<string,vector<int>> singleLineMapSeedKToIdx;
        // we also create an empty set of unique pre "vetted" Kmers in this line.
        set<string> uniqueKmersInLine;
        // populating smap for this line.
        findSeedPattern(line, singleLineMapSeedKToIdx, seedK);
        // for every smer in the map, if it appears more then once in the line:
        // we check if it can be expanded to a kmer and if so that kmer is added to the set of unique Kmers
        for (const auto& [smer, idxs] : singleLineMapSeedKToIdx){
            if (singleLineMapSeedKToIdx[smer].size() > 1){
                expandSeedToKmer(line, smer, idxs, minK, uniqueKmersInLine, legitimateSpacer, strict, maxK);
            }
        }
        // for statistics 
        if (!uniqueKmersInLine.empty()){
            numReadsWithRepeats++;
        }
        // for every unique kmer found on the line:
        // 1) we Determine the canonical form of the K-mer to ensure consistency in representation 
        //     * choosing between reverse complement and the kmer to not add both to the global kmer map
        // 2) then we increment the count of the one that was selected in the global map (or adding it for the first time) 
        for (const auto& uniqueKmer : uniqueKmersInLine) {
            string uniqueKmerOrReverse = pickKey(uniqueKmer);
            globalKmerMap[uniqueKmerOrReverse]++;
        }
    }
    // check for 0 division (process terminated before it started)
    if (progressCounter == 0){
        logFile << "no lines were processed" << endl;
        cerr << "no lines were processed" << endl;
       return;
    }
    // calculate stats:
    double precentReadsWithRepeat = (static_cast<double>(numReadsWithRepeats) / static_cast<double>(progressCounter)) * 100;
    stats["number_of_reads_in_file: "] = progressCounter;
    stats["number_of_reads_in_file_with_repeat: "] = numReadsWithRepeats;
    stats["precent_reads_in_file_with_repeat: "] = precentReadsWithRepeat;
}

bool skipThisLine(const string& read, double iligitimateRatio){
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
    if (ratio >= iligitimateRatio) {
        // Skip expansions for this line
        return true; 
    }
    return false;
}

// this function populates the smap
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

// this function populates the unique kmer set
void expandSeedToKmer(const string& line, string smer, vector<int> smerIdxVect , int minK, set<string>& uniqueKmersInLine,int legitimateSpacer, bool strict, int maxK){
    // Track unique K-mers for this particular smer on this particular line
    set<string> UniqueKmersFromSmer; 
    // we iterate over the index vector of the appearences of this smer, each index gets a turn to be the potential kmer.
    for (int potentialKmer = 0; potentialKmer < smerIdxVect.size(); potentialKmer++){
        // we initilize an empty "best fitting kmer for this index of choice"
        string bestK = "";
        // then we iterate over all the other indecies comparing the smer appearance there to our "Potential kmer"
        for (int comparisionKmer = 0; comparisionKmer < smerIdxVect.size(); comparisionKmer++ ){
            // we verify that we are no trying to expand the same location
            if (comparisionKmer == potentialKmer){
                continue;
            }
            // we set the indecies of the start and end at our potential kmer and our comparision
            int idxStartPotential = smerIdxVect.at(potentialKmer);
            int idxEndPotential = smerIdxVect.at(potentialKmer) + smer.length() - 1;
            
            int idxStartCompare = smerIdxVect.at(comparisionKmer);
            int idxEndCompare = smerIdxVect.at(comparisionKmer) + smer.length() - 1;
            
            // we set flags = true to represent are the 2 occurances equal at the start and end.
            bool equalAtStart = true;
            bool equalAtEnd = true;
            // the current Kmer is a size of an smer
            int kmerLen = smer.length();
            // while the 2 locations are equal at the start or end and the indecies are not overlapping:
            while((equalAtStart || equalAtEnd) && 
                   notOverlapping(idxStartPotential,idxStartCompare,idxEndPotential,idxEndCompare,legitimateSpacer) &&
                   kmerLen < maxK){
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
                kmerLen = idxEndPotential - idxStartPotential + 1;    
            } 
            // if they reached a point where either Kmers overlap, or hit any of the bounds, throw them out
            // if we are strict we also throw out the entire line
            if (!notOverlapping(idxStartPotential,idxStartCompare,idxEndPotential,idxEndCompare,legitimateSpacer) || 
                idxEndPotential >= (line.length() - 1) || 
                idxEndCompare >= (line.length() - 1) || 
                idxStartPotential <= 0 || 
                idxStartCompare <= 0 ||
                kmerLen >= maxK){ 
                    if(strict){
                        UniqueKmersFromSmer.clear();
                        return;
                    }
                    continue; 
            }
            // now if we reached this part of the code we are sure we have reached max expansion:
            // 1) the nucleotids dont match anymore
            // 2) they haven't reached a point of hitting boundaries or overlapping  (overlap defined as include spacer)

            // generate Kmer using substring
            //if reaches here debugg statement
            // cout << "accepted Kmer" << endl;
            string kmer = line.substr(idxStartPotential, kmerLen);
            
            // if this potential kmer follows requirements :
            // (is indeed a kmer and is longer then the best kmer so far)
            // we make it the best kmer
            if(kmerLen >= minK && kmerLen > bestK.length()){
                bestK = kmer;
            }
        }
        // once we have compared every instance of the "potential kmer" smer:
        // if we got a Best kmer we add it to the unique Kmers from smer Set
        if(bestK.size() > 0){
            UniqueKmersFromSmer.insert(bestK);
        }
            
    }
    // once we have compared every smer to all the other Smers:
    // We add the best match for each instance (a kmer) to the unique Kmers in line Set.
    for (const auto& uniqueKmer : UniqueKmersFromSmer) {
        uniqueKmersInLine.insert(uniqueKmer);
    }
}
