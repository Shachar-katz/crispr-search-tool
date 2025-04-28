#include "functions.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"

// step one functions bellow:

void findKmersInFile(MultiFormatFileReader& fileReader, 
                     unordered_map<string,int>& globalKmerMap, 
                     int seedK, 
                     int minK, 
                     int minLegitimateSpacer,
                     int maxLegitimateSpacer,
                     int horizon,
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
    // loop over every read
    while (fileReader.getNextLine(line)) {
        // statistics and progress managment:
        progressCounter++;
        if (progressCounter % interval == 0){
            cout << "Procession line: " << progressCounter << endl;
            logFile << "Procession line: " << progressCounter << endl;
            // also output faulty lines
            if (faultyLine > 0){
                logFile << faultyLine << " faulty lines" << endl;
                cerr << faultyLine << " faulty lines" << endl;
            }
            faultyLine = 0;
        }
        // check for faulty lines
        if ((preStrict == 1 && skipThisLine(line, 0.5)) || line.length() <= (2 * minK + minLegitimateSpacer + 2)){
            faultyLine++;
            continue;
        }
        
        // every line we create and populate an smap that maps from an smer to vect of indecies in the line.
        unordered_map<string,vector<int>> singleLineMapSeedKToIdx;
        findSeedPattern(line, singleLineMapSeedKToIdx, seedK);

        // we also create an empty map of unique repeats -> set of their positions in this line.
        unordered_map<string,set<int>> uniqueKmersInLine;

        // iterate over line, at each point: generate an Smer -> retrive it's index vector 
        // we check if it can be expanded to a kmer and if so that kmer is added to the set of unique Kmers
        for (int i = 0; i < line.length(); i++){
            string smer = line.substr(i,seedK);
            auto& idxs = singleLineMapSeedKToIdx[smer];
            if (idxs.size() > 1){
                expandSeedToKmer(line, smer, i, idxs, minK, uniqueKmersInLine, minLegitimateSpacer, maxLegitimateSpacer, strict, maxK, horizon);
            }
        }
        // for statistics 
        if (!uniqueKmersInLine.empty()){
            numReadsWithRepeats++;
        }
        // for every unique kmer found on the line:
        // 1) we Determine the canonic form of the K-mer to ensure consistency in representation 
        //     * choosing between reverse complement and the kmer to not add both to the global kmer map
        // 2) then we 
        for (const auto& [uniqueKmer, positions] : uniqueKmersInLine) {
            string cannonizedKmer = pickKey(uniqueKmer);
            globalKmerMap[cannonizedKmer] += positions.size();
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
bool notOverlapping(int startIdx, int idxStartCompare, int idxEnd, int idxEndCompare, int minLegitimateSpacer, int maxLegitimateSpacer)
{
    int spacing;
    if (idxEnd < idxStartCompare) {
        spacing = idxStartCompare - idxEnd -1;
        
    } 
    else if (idxEndCompare < startIdx) {
        spacing = startIdx - idxEndCompare - 1;
    }
    else {
        return false;
    }
    return (spacing >= minLegitimateSpacer && spacing <= maxLegitimateSpacer);
}

// this function populates the unique kmer set
void expandSeedToKmer(const string& line, 
                      string smer, 
                      int startIdx, 
                      vector<int> smerIdxVect, 
                      int minK, 
                      unordered_map<string, set<int>>& uniqueKmersInLine,
                      int minLegitimateSpacer, 
                      int maxLegitimateSpacer, 
                      bool strict, 
                      int maxK, 
                      int horizion){
    // Track unique K-mers for this particular smer on this particular line
    // create range
    int lowerBound = startIdx - horizion;
    int upperBound = startIdx + horizion;

    auto lowerIt = lower_bound(smerIdxVect.begin(), smerIdxVect.end(), lowerBound);
    auto upperIt = upper_bound(smerIdxVect.begin(), smerIdxVect.end(), upperBound);

    pair<string,int> bestK = {"", -1};

    // then we iterate over all the other indecies comparing the smer appearance there to our "Potential kmer"
    for(auto it = lowerIt; it != upperIt; it++){
        // cout << "entered loop" << endl; // debugg
        int comparisionPos = distance(smerIdxVect.begin(), it);
        if (comparisionPos == startIdx){ continue; }

        // we set the indecies of the start and end at our potential kmer and our comparision
        int idxEnd = startIdx + smer.length() - 1;
        
        int idxStartCompare = smerIdxVect.at(comparisionPos);
        int idxEndCompare = smerIdxVect.at(comparisionPos) + smer.length() - 1;
        
        // we set flags = true to represent are the 2 occurances equal at the start and end.
        bool equalAtStart = true;
        bool equalAtEnd = true;
        // the current Kmer is a size of an smer
        int kmerLen = smer.length();
        // while the 2 locations are equal at the start or end and the indecies are not overlapping:
        while((equalAtStart || equalAtEnd) && 
                notOverlapping(startIdx, idxStartCompare, idxEnd, idxEndCompare, minLegitimateSpacer, maxLegitimateSpacer) &&
                kmerLen < maxK){
            // as long as the indecies would valid at the start if we decremented and they are still equal at the start:
            if (equalAtStart && startIdx > 0 && idxStartCompare > 0){
                // if the "next" nucleotid is equal on both occurances we decrement the position and continue expending
                if (line[startIdx - 1] == line[idxStartCompare - 1]){
                    startIdx-- ;
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
            if(!notOverlapping(startIdx,idxStartCompare,idxEnd,idxEndCompare,minLegitimateSpacer, maxLegitimateSpacer)) {break;}
            // as long as the indecies would be valid at the end if we incremented and they are still equal at the end:
            if (equalAtEnd && (idxEnd + 1) < line.length() && (idxEndCompare + 1) < line.length()){
                // if the "next" nucleotid is equal on both occurances we increment the position and continue expending
                if (line[idxEnd + 1] == line[idxEndCompare + 1]){
                    idxEnd++ ;
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
            kmerLen = idxEnd - startIdx + 1;    
        } 
        // if they reached a point where either Kmers overlap, or hit any of the bounds, throw them out
        // if we are strict we also throw out the entire line
        if (!notOverlapping(startIdx,idxStartCompare,idxEnd,idxEndCompare,minLegitimateSpacer, maxLegitimateSpacer) || 
            idxEnd >= (line.length() - 1) || 
            idxEndCompare >= (line.length() - 1) || 
            startIdx <= 0 || 
            idxStartCompare <= 0 ||
            kmerLen >= maxK){ 
                if(strict){
                    return;
                }
                continue; 
        }
        // now if we reached this part of the code we are sure we have reached max expansion:
        // 1) the nucleotids dont match anymore
        // 2) they haven't reached a point of hitting boundaries or overlapping  (overlap defined as include spacer)
        
        // if this potential kmer follows requirements :
        // (is indeed a kmer and is longer then the best kmer so far)
        // we make it the best kmer
        if(kmerLen >= minK && kmerLen > bestK.first.length()){
            string kmer = line.substr(startIdx, kmerLen);
            bestK.first = kmer;
            bestK.second = startIdx;
        }   
    }
    // once we have compared every smer to all the other Smers:
    // We add the best match for each instance (a kmer) to the unique Kmers in line Set.
    if(bestK.second != -1){
        auto& positions = uniqueKmersInLine[bestK.first];
        positions.insert(bestK.second);
    }
}
