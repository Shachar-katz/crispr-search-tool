//
//  functions.hpp
//  repeatesSearchProject
//
//  Created by Sarah Katz on 10/23/24.
//

#ifndef functions_hpp
#define functions_hpp

#include <stdio.h>
#include <cctype>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <algorithm>
#include "DynamicBinsClass.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"
#include "arrayClass.h"

using namespace std;


// step one functions:
bool isInputFileValid(ifstream& inFS, 
                      string fileName);

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
                     int maxK);

bool skipThisLine(const string& read, double iligitimateRatio);

void findSeedPattern(string line, unordered_map<string,vector<int>>& singleLineMapSeedKToIdx, int seedK);

bool notOverlapping(int idxStartPotential, 
                    int idxStartCompare, 
                    int idxEndPotential, 
                    int idxEndCompare, 
                    int minLegitimateSpacer);

int expandSeedToKmer(const string& line, 
                      string smer, 
                      int startIdx, 
                      vector<int> smerIdxVect, 
                      int minK, 
                      unordered_map<int, string>& PosToKmersInLine,
                      int minLegitimateSpacer, 
                      bool strict, 
                      int maxK, 
                      int horizion);

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
                         unordered_map<int,string>& posToKmerInLine);

void smoothRepScore(const vector<double>& inLineRepetitionScores, 
                    int smoothingWindow,
                    vector<double>& smoothedScores);

void createExclusionMask(const vector<double>& smoothedScores,
                         int exclusionWindow,
                         vector<bool>& excludedSegments);

bool excludeLine(const vector<double>& smoothedScores,
                int exclusionMinWindows); // debugg


//step two functions:

void catalogToSAndKMaps(ifstream& inCatalog, unordered_map<string,vector<string>>& smap, 
                        unordered_map<string,data_t>& kmap, int seedK, ofstream& logFile);
void findSmersVect(string shortestKmer, vector<string>& smerVect, int seedK);

void makePotentialRelationsSet(const unordered_map<string,vector<string>>& smap, 
                               set<pair<string,string>>& potentialRelationSet, ofstream& logFile);
void chooseShortestK (const pair<string, string>& kPair, string& shortestK, string& otherK);

void binRelatives(const string& shortestK, const string& otherK, DynamicBins& bins);

void verifyRelation(const unordered_map<string,vector<string>>& smap, const set<pair<string,string>>& potentialRelationSet, 
                    int seedK, DynamicBins& bins, ofstream& logFile, int alpha = 1);

string kmerCompetition(const unordered_map<string,data_t>& kmap, string currentRep, string auditioningKmer, ofstream& logFile);

string tieBreaker(string currentRep, string auditioningKmer);

void binSingles(const unordered_map<string,vector<string>>& smap, DynamicBins& bins, ofstream& logFile);

void selectReps(unordered_map<int, string>& provisionalRepList, const unordered_map<int,
                vector<string>>& reverseBins, const unordered_map<string,data_t>& kmap, ofstream& logFile);

void selectRepsWeight(unordered_map<int, string>& provisionalRepList, 
                      const unordered_map<int, 
                      vector<string>>& reverseBins, 
                      const unordered_map<string,data_t>& kmap, 
                      int seedK, 
                      ofstream& logFile);
                    
inline string findRepUsingWeight(const vector<string>& kVect, 
                                 const unordered_map<string,data_t>& kmap, 
                                 int seedK);

inline int distance(string kmerI, string kmerJ , int seedK);

unordered_map<int, string> reCannonization(const unordered_map<int, string>& provisionalRepList, const DynamicBins& bins, ofstream& logFile);

void creatingOutputMap(unordered_map<string,data_t>& outputMap, unordered_map<string,data_t>& binsOutputMap, 
                       const unordered_map<int,string>& choosenReps, const unordered_map<int,
                       vector<string>>& reverseBins, const unordered_map<string,data_t>& kmap);

int palindromicScore(string kmer, int alpha = 0);

void validateBins(const unordered_map<int, string>& provisionalRepList, const DynamicBins& bins, 
                  const unordered_map<int,vector<string>>& reverseBins, ofstream& logFile);

inline string largeClusterPartition(const unordered_map<string,data_t>& kmap, const vector<string>& kVect, int seedK);

// bool isKmer(string maybeKmer);

// step three functions:
int buildSmap(ifstream& inCatalog, unordered_map<string,Kmap_t>& smap, int seedK, int minPalindromic = 0);

void findKmersInFileWithSmap(MultiFormatFileReader& fileReader, unordered_map<string,data_t>& globalKmerMap, 
                            unordered_map<string,Kmap_t>& smap, int seedK, unordered_map<string,double>& stats, ofstream& logFile, int minLegitimateSpacer, int minK, int interval);
void expandSeedToKmerWithSmap(const string& line, const string& smer, int& idxInLine, unordered_map<string,data_t>& globalKmerMap, unordered_map<string,Kmap_t>& smap, bool& activeLine, unordered_map<string,int>& kmerToIdxInLine);

bool willSelfOverlap(const unordered_map<string,int>& kmerToIdxInLine,int startIdexOfKmerInLine, string kmerInLine);

bool valideHeader(string header);

// Array Dump functions:
int buildKmap(ifstream& inCatalog, unordered_map<string,string>& kmerToId, int minPalindromic = 0);

void arrayIdentifior(MultiFormatFileReader& fileReader, 
                     unordered_map<string,Array>& globalArrayVect, 
                     unordered_map<string,Kmap_t>& smap,
                     unordered_map<string,string>& kmerToId, 
                     int seedK, 
                     unordered_map<string,double>& stats, 
                     ofstream& logFile, 
                     int minLegitimateSpacer, 
                     int maxLegitimateSpacer, 
                     int minK, 
                     int interval,
                     int maxMismatches);

string expandSeedToKmer(const string& line, 
                        const string& smer, 
                        int& idxInLine, 
                        unordered_map<string,Kmap_t>& smap, 
                        bool& activeLine, 
                        int& tempStartIdx, 
                        int maxMismatches, 
                        int& numMissmatches);

inline void findSmerSet(string kmer, unordered_set<string>& smerSet, int seedK){
    for (int j = 0; j <= (kmer.size() - seedK); j++){
        string smerForSearch = kmer.substr(j,seedK);
        smerSet.insert(smerForSearch);
    }
}

inline bool isKmerMatch(const string& line,  
                        int start,
                        int end, 
                        const unordered_set<string>& smerSet, 
                        int seedK, 
                        int& missmatches, 
                        int maxMismatches = 0);

void constructRepeatMap(const unordered_map<string,Array>& globalArrayMap, 
                        unordered_map<string,RepeatData>& repeatMap);


#endif /* functions_hpp */
