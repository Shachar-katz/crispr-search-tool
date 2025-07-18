#include "functions.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"

// step one functions bellow:

void findKmersInFile(MultiFormatFileReader& fileReader, 
                     unordered_map<string,int>& globalKmerMap, 
                     int seedK, 
                     int minK, 
                     int minLegitimateSpacer,
                     int horizon,
                     int segmentSize,
                     int smoothingWindow,
                     unordered_map<string,double>& stats, 
                     bool strict, 
                     bool preStrict, 
                     ofstream& logFile,
                     int interval,
                     int maxK,
                     const string& currentIdentifier,
                     unordered_map<string,set<string>>& kmerToIdentifiers)
    {
    // line variable temporerally holds the reads
    string line;
    // statistics Vars
    int progressCounter = 0;
    int numReadsWithRepeats = 0;
    int faultyLine = 0;
    ofstream dump;
    string dumpName = "/Users/sarahkatz/relman_lab/L1_V5/step_1_repeat_dump";
    dump.open(dumpName);
    dump << "repeat" << '\t' << "line_num" <<  '\t' << "pos" << endl;
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

        // we map the entire line to segments of equal size that each recive a repetition score defined mathmatically 
        vector<double> inLineRepetitionScores;
        unordered_map<int, string> posToKmersInLine;
        generateRepeatition(line, segmentSize, seedK, minK, maxK, minLegitimateSpacer, horizon, smoothingWindow, strict, singleLineMapSeedKToIdx, inLineRepetitionScores, posToKmersInLine);
        vector<double> smoothScores;
        smoothRepScore(inLineRepetitionScores, smoothingWindow, smoothScores); // do i need this?
        for (auto& [pos, repeat]: posToKmersInLine){
            dump << repeat << '\t' << progressCounter <<  '\t' << pos << endl;
        }

        // if (excludeLine(smoothScores, 30) == true) { continue; } // debugg
        // then we mask the segments with a repition score of over 0.6 and their neighboring 2 segments on each side
        vector<bool> maskedSegments;
        int exclusionWindow = 3; // not needed i think throwing out the line should do !!
        createExclusionMask(smoothScores, exclusionWindow, maskedSegments); //  a change!!!
        // then we populate a map of unique repeats -> set of their positions in this line based only on repeats found in the unmasked areas
        unordered_map<string,set<int>> uniqueKmersInLine;
        for(int i = 0; i < maskedSegments.size(); i++){
            if (maskedSegments[i] == true) { continue; }
            int segment = i * segmentSize;
            int segmentEnd = segment + segmentSize;
            for (int j = segment; j < segmentEnd; j++){
                if (posToKmersInLine.count(j) == 0) { continue; }
                string kmer = posToKmersInLine.at(j);
                auto& positions = uniqueKmersInLine[kmer];
                positions.insert(j);
            }
        }
        // for (auto& [pos, repeat]: posToKmersInLine){
        //     auto& positions = uniqueKmersInLine[repeat];
        //     positions.insert(pos);
        // }

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
            kmerToIdentifiers[cannonizedKmer].insert(currentIdentifier); // for tracking which were found where
        }
    }
    dump.close();
    // check for 0 division (process terminated before it started)
    if (progressCounter == 0){
        logFile << "no lines were processed" << endl;
        cerr << "no lines were processed" << endl;
       return;
    }
    // calculate stats:
    stats["number_of_reads_in_file: "] += progressCounter;
    stats["number_of_reads_in_file_with_repeat: "] += numReadsWithRepeats;
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
    for (int j = 0; j <= (line.size() - seedK); j++){
        string key = line.substr(j,seedK);
        if (key.size() < seedK){
            break; // throw error
        }
        singleLineMapSeedKToIdx[key].emplace_back(j);
    }
}

// this function recives start and end location of 2 Kmers in expansion and 
// varifies that they are not overlapping
bool notOverlapping(int startIdx, int idxStartCompare, int idxEnd, int idxEndCompare, int minLegitimateSpacer)
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
    return (spacing >= minLegitimateSpacer);
}

// this function populates the unique kmer set
int expandSeedToKmer(const string& line, 
                      string smer, 
                      int startIdx, 
                      vector<int> smerIdxVect, 
                      int minK, 
                      unordered_map<int, string>& PosToKmersInLine,
                      int minLegitimateSpacer, 
                      bool strict, 
                      int maxK, 
                      int horizion){
    // Track unique K-mers for this particular smer on this particular line
    // create range
    int lowerBound = startIdx - horizion;
    int upperBound = startIdx + horizion;

    auto lowerIt = lower_bound(smerIdxVect.begin(), smerIdxVect.end(), lowerBound);
    auto upperIt = upper_bound(smerIdxVect.begin(), smerIdxVect.end(), upperBound);

    string bestK = "";
    int bestKPos = -1;

    // then we iterate over all the other indecies comparing the smer appearance there to our "Potential kmer"
    for(auto it = lowerIt; it != upperIt; it++){
        if (*it == startIdx){ continue; }

        // we set the indecies of the start and end at our potential kmer and our comparision
        int startIdxCopy = startIdx;
        int idxEnd = startIdx + smer.length() - 1;
        
        int idxStartCompare = *it;
        int idxEndCompare = idxStartCompare + smer.length() - 1;
        
        // we set flags = true to represent are the 2 occurances equal at the start and end.
        bool equalAtStart = true;
        bool equalAtEnd = true;
        // the current Kmer is a size of an smer
        int kmerLen = smer.length();
        // while the 2 locations are equal at the start or end and the indecies are not overlapping:
        while((equalAtStart || equalAtEnd) && 
                notOverlapping(startIdxCopy, idxStartCompare, idxEnd, idxEndCompare, minLegitimateSpacer) &&
                kmerLen < maxK){
            // as long as the indecies would valid at the start if we decremented and they are still equal at the start:
            if (equalAtStart && startIdxCopy > 0 && idxStartCompare > 0){
                // if the "next" nucleotid is equal on both occurances we decrement the position and continue expending
                if (line[startIdxCopy - 1] == line[idxStartCompare - 1]){
                    startIdxCopy-- ;
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
            if(!notOverlapping(startIdxCopy,idxStartCompare,idxEnd,idxEndCompare,minLegitimateSpacer)) {break;}
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
            kmerLen = idxEnd - startIdxCopy + 1;    
        } 
        // if they reached a point where either Kmers overlap, or hit any of the bounds, throw them out
        // if we are strict we also throw out the entire line
        if (!notOverlapping(startIdxCopy,idxStartCompare,idxEnd,idxEndCompare,minLegitimateSpacer) || 
            idxEnd >= (line.length() - 1) || 
            idxEndCompare >= (line.length() - 1) || 
            startIdxCopy <= 0 || 
            idxStartCompare <= 0 ||
            kmerLen >= maxK){ 
                if(strict){
                    return 0;
                }
                continue; 
        }
        // now if we reached this part of the code we are sure we have reached max expansion:
        // 1) the nucleotids dont match anymore
        // 2) they haven't reached a point of hitting boundaries or overlapping  (overlap defined as include spacer)
        
        // if this potential kmer follows requirements :
        // (is indeed a kmer and is longer then the best kmer so far)
        // we make it the best kmer
        if(kmerLen >= minK && kmerLen > bestK.length()){
            string kmer = line.substr(startIdxCopy, kmerLen);
            bestK = kmer;
            bestKPos = startIdxCopy;
        }   
    }
    // once we have compared every smer to all the other Smers:
    // We add the best match for each instance (a kmer) to the unique Kmers in line Set.
    if(bestKPos != -1){
        auto& repeat = PosToKmersInLine[bestKPos];
        repeat = bestK;
        return 1;
    }
    return 0;
}

// calculates repeatition score that is defined as: the number of positions from which an smer successfully extended into a repeat
// returns a smoothed rep score ( bins and calculates a weighted avg of the 2 bins before and afer a bin + the bin itself )
void generateRepeatition(const string& line,
                         int segmentSize,
                         int seedK, 
                         int minK,
                         int maxK,
                         int minLegitimateSpacer,
                         int horizon,
                         int smoothingWindow,
                         bool strict,
                         const unordered_map<string,vector<int>>& singleLineMapSeedKToIdx, 
                         vector<double>& inLineRepetitionScores,
                         unordered_map<int,string>& posToKmerInLine)
{

    int currSegment = 0;
    double currRepitionScore = 0;
    for (int i = 0; i <= (line.length() - seedK) ; i++){
        if (i % segmentSize == 0 && i != 0){
            inLineRepetitionScores.emplace_back(currRepitionScore / segmentSize);
            currSegment = i; 
            currRepitionScore = 0.0;
        }
        string smer = line.substr(i,seedK);
        auto& idxs = singleLineMapSeedKToIdx.at(smer);
        if (idxs.size() > 1){
            currRepitionScore += expandSeedToKmer(line, smer, i, idxs, minK, posToKmerInLine, minLegitimateSpacer, strict, maxK, horizon);
        }
    }
    if (currSegment < line.length()) {
        int remainingPositions = (line.length() - seedK + 1) - currSegment; // Positions actually processed
        double finalScore = (remainingPositions > 0) ? currRepitionScore / remainingPositions : 0.0;
        inLineRepetitionScores.emplace_back(finalScore);
    }

}

void smoothRepScore(const vector<double>& inLineRepetitionScores, 
                    int smoothingWindow,
                    vector<double>& smoothedScores){
    
    int numSegments = static_cast<int>(inLineRepetitionScores.size());
    for (int j = 0; j < numSegments; j++)
    {
        int lowIdx = std::max(j - smoothingWindow, 0);
        int highIdx = std::min(j + smoothingWindow, numSegments - 1);
        double smoothedScore = 0.0;
        double sumA = 0.0;
        for (int i = lowIdx; i <= highIdx; i++){
            int aI = (2 * smoothingWindow + 1) - abs(i - j);
            double weightI = static_cast<double>(aI);
            smoothedScore += weightI * inLineRepetitionScores[i];
            sumA += aI;
        }
        smoothedScores.emplace_back( smoothedScore / sumA );
    }
}

void createExclusionMask(const vector<double>& inLineRepetitionScores,
                         int exclusionWindow,
                         vector<bool>& excludedSegments) {
    int numSegments = inLineRepetitionScores.size();
    excludedSegments.resize(numSegments, false);
    // Find and exclude segments adjacent to high repetition areas (> 0.6)
    int i = 0;
    while (i < numSegments) {
        if (inLineRepetitionScores[i] < 0.55)
        { 
            i++;
            continue; 
        }
        int lowIdx = std::max(i - exclusionWindow, 0);
        int highIdx = std::min(i + exclusionWindow, numSegments - 1);
        for (int j = lowIdx; j <= highIdx; j++){
            excludedSegments[j] = true;
        }
        i += exclusionWindow;
    }
}

bool excludeLine(const vector<double>& smoothedScores,
                int exclusionMinWindows) {
    int numSegments = smoothedScores.size();
    // Find and exclude segments adjacent to high repetition areas (> 0.6)
    int overLimCount = 0;
    for (int i = 0; i < numSegments; i++)
    {
        if (smoothedScores[i] > 0.6)
        { 
            overLimCount++;
        }
        if (overLimCount >= exclusionMinWindows){
            return true;
        }
    }
    return false;
}