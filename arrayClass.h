#pragma once
#include <iostream>
#include <string>
#include <vector>
#include "inlines.h"
using namespace std;


// class Array{
// private:
//     string repeat;
//     string repeatId;
//     vector<int> inLineCoordinatesVect;
//     vector<string> array;
//     vector <int> repeatIdxToNumMissmatches;
//     // add spacer object?
//     int arrayLen = 0; // might not need?
//     int spacerLen = 0;
//     int numSpacers;
//     bool invalid = false;
// public:
//     Array(){}
//     Array(string repeat)
//     {
//         this->repeat = repeat;
//     }
//     ~Array(){}
//     void openArray(string repeat, int startIdxInLine, string repeatId, int numMissmatches)
//     {
//         this->repeat = repeat;
//         this->repeatId = repeatId;
//         addRepeat(startIdxInLine, numMissmatches);
//         if (inLineCoordinatesVect.size() == 2){ // set spacer length
//             this->spacerLen = inLineCoordinatesVect[1] - inLineCoordinatesVect[0]; 
//         }
//     }
//     void addRepeat(int startIdxInLine, int numMissmatches)
//     {   
//         int endIdxInLine = startIdxInLine + repeat.length();
//         inLineCoordinatesVect.push_back(startIdxInLine);
//         inLineCoordinatesVect.push_back(endIdxInLine);
//         repeatIdxToNumMissmatches.push_back(numMissmatches);
//     }
//     bool stillValid(int startIdxInLine, int maxValidSpacer, int minValidSpacer, bool perfectMatch){
//         int spacerLen = startIdxInLine - inLineCoordinatesVect.back();
//         if ((spacerLen <= maxValidSpacer && (startIdxInLine > (inLineCoordinatesVect.back() + minValidSpacer) || perfectMatch))){ return true; }
//         return false;
//     }
//     // void setInvalid(bool invalid) { hitInvalid; }
//     string getRepeat() const { return repeat; }
//     string getRepeatId() const { return repeatId; }
//     int getNumSpacers() const { return numSpacers; }
//     int getArrayLen() const { return arrayLen; }
//     int getSpacerLen() const { return spacerLen; }
//     int getStartPos() const { return inLineCoordinatesVect[0]; }
//     int getEndPos() const { return (inLineCoordinatesVect[0] + arrayLen); }
//     int getStartError() const { return inLineCoordinatesVect[inLineCoordinatesVect.size() - 3]; }
//     bool hitInvalid() const { return invalid; }
//     vector<string> getArrayVect() const { return this->array; }
//     string getArrayStr() 
//     {
//         string array;
//         for (int i = 0; i < this->array.size(); i++)
//         {
//             array += this->array[i] + "|";
//         }
//         return array;
//     }
//     vector<int> getErrorVect() const { return repeatIdxToNumMissmatches; }
//     bool closeArray(string line){
//         if (inLineCoordinatesVect.size() < 4) { return false; }
//         arrayLen = inLineCoordinatesVect.back() - inLineCoordinatesVect[0];
//         for (int i = 0; i < inLineCoordinatesVect.size() - 1; i++)
//         {
//             int segmentLen = inLineCoordinatesVect[i + 1] - inLineCoordinatesVect[i];
//             string arraySegment = line.substr(inLineCoordinatesVect.at(i), segmentLen);
//             array.push_back(arraySegment);
//             // potentially populate repeat in array coordinate vect
//         }
//         // calculate num spacers
//         int numRepeats = inLineCoordinatesVect.size() / 2;
//         numSpacers = numRepeats - 1;
//         return true;
//     }


// };

// class LineArrayHandler{
// private:
//     string line;
//     int maxAllowedSpacer;
//     int minAllowedSpacer;
//     vector<Array> lineArrayVect;
//     bool activeArray = false;
//     Array tempArray;
//     int startErrorField;
// public:
//     LineArrayHandler(const string& line, int maxAllowedSpacer, int minAllowedSpacer)
//     {
//         this->line = line;
//         this->maxAllowedSpacer = maxAllowedSpacer;
//         this->minAllowedSpacer = minAllowedSpacer;
//     }
//     ~LineArrayHandler(){}
//     void manageState(string repeat, int startIdxInLine, string repeatId, int numMissmatches)
//     {
//         if(!activeArray){
//             tempArray.openArray(repeat, startIdxInLine, repeatId, numMissmatches);
//             activeArray = true;
//         }
//         // else if(tempArray.getRepeat() == repeat){
//         else if(areRepeatsTheSame(tempArray.getRepeat(), repeat, 10, 15)){
//             bool perfectMatch = areRepeatsTheSame(tempArray.getRepeat(), repeat, 10, 3) && abs(static_cast<int>(tempArray.getRepeat().length() - repeat.length())) < 4;
//             bool wasExpanded = this->expandArray(startIdxInLine, numMissmatches, perfectMatch);
//             if (!wasExpanded){
//                 tempArray.openArray(repeat, startIdxInLine, repeatId, numMissmatches);
//                 activeArray = true;
//             }
//         }
//         // else if(tempArray.getRepeat() != repeat){
//         else{
//             this->uploadArray(); // should I use return value for anything?
//             tempArray.openArray(repeat, startIdxInLine, repeatId, numMissmatches);
//             activeArray = true;
//         }
//     }
//     bool expandArray(int startIdxInLine, int numMissmatches, bool perfectMatch)
//     {
//         if (tempArray.stillValid(startIdxInLine, maxAllowedSpacer, minAllowedSpacer, perfectMatch)){
//             tempArray.addRepeat(startIdxInLine, numMissmatches);
//             return true;
//         }
//         // if (!tempArray.hitInvalid()){
//         //     startErrorField = tempArray.getStartError();
//         //     tempArray.setInvalid(true);
//         //     return true;
//         // }
//         bool uploaded = uploadArray();
//         return false;
//     }
//     bool uploadArray()
//     {
//         bool validArray = tempArray.closeArray(line);
//         if (validArray){ 
//             lineArrayVect.push_back(tempArray);
//         }
//         tempArray = Array();
//         activeArray = false;
//         return validArray;
//     }
//     bool isActive() { return activeArray; }
//     int getSpacerLen() const { return tempArray.getSpacerLen(); }
//     vector<Array> getLineArrayVect(){ return lineArrayVect; }
//     string getCurrRepeat() const { return tempArray.getRepeat(); }
//     bool noArrays() {
//         if (lineArrayVect.empty()){
//             return true;
//         }
//         return false;
//     }
    
// };

class Array{
private:
    string repeat;
    string repeatId;
    vector<int> inLineCoordinatesVect;
    vector<int> inArrayRepeatCoordinatesVect; // might not need?
    vector<string> array;
    vector<int> repeatIdxToNumMissmatches;
    int arrayLen = 0;
    int numSpacers;
    int avgSpacerLen = 0;  // NEW: tracks average spacer length
    
    int calculateAvgSpacerLen()
    {
        if (inLineCoordinatesVect.size() < 4){  // Need at least 2 repeats
            return 0;
        }
        
        int totalSpacerLen = 0;
        int numSpacers = 0;
        
        // Iterate over spacers (between repeat pairs)
        // Coords: [R1_start, R1_end, R2_start, R2_end, ...]
        // Spacer is from R1_end to R2_start
        for (size_t i = 1; i < inLineCoordinatesVect.size() - 1; i += 2){
            int spacerStart = inLineCoordinatesVect[i];
            int spacerEnd = inLineCoordinatesVect[i + 1];
            int spacerLen = spacerEnd - spacerStart;
            totalSpacerLen += spacerLen;
            numSpacers++;
        }
        
        if (numSpacers == 0) return 0;
        return totalSpacerLen / numSpacers;
    }

public:
    Array(){}
    Array(string repeat)
    {
        this->repeat = repeat;
    }
    ~Array(){}
    
    void openArray(const string& repeat, int startIdxInLine, const string& repeatId, int numMissmatches)
    {
        this->repeat = repeat;
        this->repeatId = repeatId;
        addRepeat(startIdxInLine, numMissmatches);
    }
    
    void addRepeat(int startIdxInLine, int numMissmatches)
    {   
        int endIdxInLine = startIdxInLine + repeat.length();
        inLineCoordinatesVect.push_back(startIdxInLine);
        inLineCoordinatesVect.push_back(endIdxInLine);
        repeatIdxToNumMissmatches.push_back(numMissmatches);
        
        // UPDATE AVERAGE SPACER LENGTH
        avgSpacerLen = calculateAvgSpacerLen();
    }
    
    bool stillValid(int startIdxInLine, int maxValidSpacer, int minValidSpacer){
        int spacerLen = startIdxInLine - inLineCoordinatesVect.back();
        if (spacerLen <= maxValidSpacer && startIdxInLine > (inLineCoordinatesVect.back() + minValidSpacer)){ 
            return true; 
        }
        return false;
    }
    
    string getRepeat() const { return repeat; }
    string getRepeatId() const { return repeatId; }
    int getNumSpacers() const { return numSpacers; }
    int getArrayLen() const { return arrayLen; }
    vector<string> getArrayVect() const { return this->array; }
    
    string getArrayStr() 
    {
        string arrayStr;
        for (size_t i = 0; i < this->array.size(); i++)
        {
            arrayStr += this->array[i] + "|";
        }
        return arrayStr;
    }
    
    vector<int> getErrorVect() const { return repeatIdxToNumMissmatches; }
    
    int getLastRepeatEnd() const
    {
        if (inLineCoordinatesVect.empty()){
            return -1;
        }
        return inLineCoordinatesVect.back();
    }
    
    int getAvgSpacerLen() const
    {
        return avgSpacerLen;
    }
    
    int getStartPos() const
    {
        if (inLineCoordinatesVect.empty()){
            return -1;
        }
        return inLineCoordinatesVect[0];
    }
    
    int getEndPos() const
    {
        return getLastRepeatEnd();
    }
    
    bool closeArray(const string& line){
        if (inLineCoordinatesVect.size() < 4) { 
            return false; 
        }
        
        arrayLen = inLineCoordinatesVect.back() - inLineCoordinatesVect[0];
        
        for (size_t i = 0; i < inLineCoordinatesVect.size() - 1; i++)
        {
            int segmentLen = inLineCoordinatesVect[i + 1] - inLineCoordinatesVect[i];
            string arraySegment = line.substr(inLineCoordinatesVect.at(i), segmentLen);
            array.push_back(arraySegment);
        }
        
        // calculate num spacers
        int numRepeats = inLineCoordinatesVect.size() / 2;
        numSpacers = numRepeats - 1;
        return true;
    }
};

class LineArrayHandler{
private:
    string line;
    int maxAllowedSpacer;
    int minAllowedSpacer;
    int seedK;
    int maxMismatches;
    unordered_map<int, Array> internalArrayMap;
    int internalArrayCounter;
    
    double countWeightedMismatches(const string& seq1, const string& seq2) {
        if (seq1.length() != seq2.length()){
            return (double)INT_MAX;
        }
        
        double totalMismatches = 0.0;
        int frontEnd = seedK * 3;
        int tailStart = (int)seq1.length() - (seedK * 3);
        
        if (tailStart < frontEnd) {
            frontEnd = 0;
            tailStart = seq1.length();
        }
        
        for (size_t i = 0; i < seq1.length(); i++){
            if (seq1[i] != seq2[i]){
                if ((int)i < frontEnd || (int)i >= tailStart){
                    totalMismatches += 0.1;
                }
                else{
                    totalMismatches += 1.0;
                }
            }
        }
        
        return totalMismatches;
    }
    
    void addRepeatToSmerMap(const string& repeat, 
                           unordered_map<string, unordered_map<string, vector<int>>>& smap)
    {
        for (int j = 0; j <= (int)(repeat.length() - seedK); j++){
            string smer = repeat.substr(j, seedK);
            smap[smer][repeat].push_back(j);
        }
    }

    double compareRepeatsLengthAdjusted(const string& candidate, const string& reference)
    {
        // Compare only up to the length of the shorter sequence
        int compareLength = min(candidate.length(), reference.length());
        
        double totalMismatches = 0.0;
        int frontEnd = seedK * 3;
        int tailStart = compareLength - (seedK * 3);
        
        if (tailStart < frontEnd){
            frontEnd = 0;
            tailStart = compareLength;
        }
        
        for (int i = 0; i < compareLength; i++){
            if (candidate[i] != reference[i]){
                if (i < frontEnd || i >= tailStart){
                    totalMismatches += 0.1;
                }
                else{
                    totalMismatches += 1.0;
                }
            }
        }
        
        // Penalize length difference
        int lengthDiff = abs((int)candidate.length() - (int)reference.length());
        totalMismatches += lengthDiff * 0.5; // Each bp difference costs 0.5 mismatch
        
        return totalMismatches;
    }

    double calculateWeightedSmerCoverage(const string& candidateRepeat, 
                                     const unordered_map<string, unordered_map<string, vector<int>>>& dynamicSmap,
                                     int seedK)
    {
        // Build match vector for all smers
        vector<bool> smerMatches;
        for (int i = 0; i <= (int)candidateRepeat.length() - seedK; i++){
            string smer = candidateRepeat.substr(i, seedK);
            smerMatches.push_back(dynamicSmap.find(smer) != dynamicSmap.end());
        }
        
        int totalSmers = smerMatches.size();
        if (totalSmers == 0) return 0.0;
        
        double effectiveMismatches = 0.0;
        int edgeSize = 3; // First 3 and last 3 smers
        int i = 0;
        
        while (i < totalSmers){
            if (smerMatches[i]){
                // Match - no penalty
                i++;
            }
            else{
                // Mismatch - apply weighted penalty
                bool isEdge = (i < edgeSize || i >= (totalSmers - edgeSize));
                
                if (isEdge){
                    // Edge region - only 1/6 penalty per mismatch
                    effectiveMismatches += (1.0 / 6.0);
                    i++;
                }
                else{
                    // Middle region - check for consecutive run (potential indel)
                    int runStart = i;
                    while (i < totalSmers && !smerMatches[i]){
                        i++;
                    }
                    int runLength = i - runStart;
                    
                    if (runLength == seedK){
                        // Exactly seedK consecutive mismatches = single indel
                        effectiveMismatches += 1.0;
                    }
                    else{
                        // Other lengths = full penalty
                        effectiveMismatches += runLength;
                    }
                }
            }
        }
        
        // Coverage = 1 - (effective mismatches / total)
        return 1.0 - (effectiveMismatches / (double)totalSmers);
    }

public:
    LineArrayHandler(const string& line, int maxAllowedSpacer, int minAllowedSpacer, int seedK, int maxMismatches)
    {
        this->line = line;
        this->maxAllowedSpacer = maxAllowedSpacer;
        this->minAllowedSpacer = minAllowedSpacer;
        this->seedK = seedK;
        this->maxMismatches = maxMismatches;
        this->internalArrayCounter = 0;
    }
    
    ~LineArrayHandler(){}

    int buildArrayFromFirstRepeat(const string& firstRepeat, int firstStartIdx, 
                               const string& firstRepeatId, int firstNumMismatches)
    {
        Array tempArray;
        tempArray.openArray(firstRepeat, firstStartIdx, firstRepeatId, firstNumMismatches);
        
        unordered_map<string, unordered_map<string, vector<int>>> dynamicSmap;
        addRepeatToSmerMap(firstRepeat, dynamicSmap);
        
        int referenceLength = firstRepeat.length();
        
        while (true){
            int consecutiveShortOrBad = 0; // Track declining quality
            const int MAX_DECLINING = 2;   // Stop if 2 consecutive bad repeats
            int minValidStartPos = tempArray.getLastRepeatEnd() + minAllowedSpacer;
            int maxValidStartPos = tempArray.getLastRepeatEnd() + maxAllowedSpacer;
            
            if (minValidStartPos > (int)(line.length() - seedK)){
                break;
            }
            
            bool foundMatch = false;
            int matchStartIdx = -1;
            double matchMismatches = 0;
            string matchRepeat = "";
            int matchLength = 0;
            
            cerr << "DEBUG: Scanning positions [" << minValidStartPos << "," << maxValidStartPos << "]" << endl;
            
            // Scan positions in order - take FIRST valid match
            for (int startPos = minValidStartPos; 
                startPos <= maxValidStartPos && startPos <= (int)(line.length() - seedK); 
                startPos++){
                
                int minLength = referenceLength * 0.8;
                int maxLength = min((int)(referenceLength * 1.1), (int)(line.length() - startPos)); // Allow up to 110% // Changed: cap at referenceLength
                
                // Try different lengths for this position
                for (int length = minLength; length <= maxLength; length++){
                    
                    string candidateRepeat = line.substr(startPos, length);
                    
                    // Calculate smer coverage
                    
                    double smerCoverage = calculateWeightedSmerCoverage(candidateRepeat, dynamicSmap, seedK);
                    
                    if (startPos >= 1620 && startPos <= 1650){
                        cerr << "  pos=" << startPos << " len=" << length 
                            << " coverage=" << smerCoverage << endl;
                    }
                    if (smerCoverage < 0.50){
                        continue;
                    }
                    
                    // Compare to all repeats in dynamic map
                    double bestMismatch = maxMismatches + 1.0;
                    
                    for (const auto& [smer, kmap] : dynamicSmap){
                        for (const auto& [refRepeat, indices] : kmap){
                            double mismatch = compareRepeatsLengthAdjusted(candidateRepeat, refRepeat);
                            
                            if (mismatch < bestMismatch){
                                bestMismatch = mismatch;
                            }
                            
                            if (mismatch <= maxMismatches){
                                break;
                            }
                        }
                        if (bestMismatch <= maxMismatches){
                            break;
                        }
                    }
                    
                    // Take FIRST match that passes threshold
                    if (bestMismatch <= maxMismatches){
                        matchStartIdx = startPos;
                        matchMismatches = bestMismatch;
                        matchRepeat = candidateRepeat;
                        matchLength = length;
                        foundMatch = true;
                        break; // Stop trying other lengths for this position
                    }
                }
                
                if (foundMatch){
                    break; // Stop scanning other positions - we found the first valid one
                }
            }
            
            if (foundMatch){
                bool isDeclining = false;
                if (matchLength < referenceLength * 0.5 || matchMismatches > maxMismatches * 0.7){
                    consecutiveShortOrBad++;
                    isDeclining = true;
                }
                else{
                    consecutiveShortOrBad = 0; // Reset counter on good match
                }
                
                // Stop array if quality has degraded
                if (consecutiveShortOrBad >= MAX_DECLINING){
                    cerr << "DEBUG: Quality degraded (short/high-mismatch repeats), stopping array" << endl;
                    break;
                }
                cerr << "DEBUG: Added repeat at pos " << matchStartIdx 
                    << " length=" << matchLength
                    << " with " << matchMismatches << " weighted mismatches";
                if (isDeclining){
                    cerr << " [DECLINING QUALITY]";
                }
                cerr << endl;
                    
                tempArray.addRepeat(matchStartIdx, matchMismatches);
                addRepeatToSmerMap(matchRepeat, dynamicSmap);
            }
            else{
                cerr << "DEBUG: No more repeats found" << endl;
                break;
            }
        }
        
        int arrayEndIdx = tempArray.getLastRepeatEnd();
        bool isValid = tempArray.closeArray(line);
        
        if (isValid){
            internalArrayMap[internalArrayCounter] = tempArray;
            internalArrayCounter++;
            cerr << "DEBUG: Stored valid array" << endl;
        }
        
        return arrayEndIdx;
    }

    
    vector<Array> extractValidArrays()
    {
        vector<Array> arrayVector;
        for (const auto& [id, array] : internalArrayMap){
            arrayVector.push_back(array);
        }
        return arrayVector;
    }
};

ostream& operator<<(ostream& out, const Array& obj);

class RepeatData{
public:
    string arrayID;
    string repeatID;
    int repeatIdx;
    int numMissmatches;
    int repeatLen;
};

ostream& operator<<(ostream& out, const RepeatData& obj);

class ArrayPositionData{
public:
    string readID;
    int startPos;
    int endPos;
    int readLen;
};

ostream& operator<<(ostream& out, const ArrayPositionData& obj);