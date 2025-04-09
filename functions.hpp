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
#include <stdexcept>
#include <algorithm>
#include "DynamicBinsClass.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"

using namespace std;


// step one functions:
bool isInputFileValid(ifstream& inFS, string fileName);

void findKmersInFile(MultiFormatFileReader& fileReader, unordered_map<string,int>& globalKmerMap, int seedK, 
                     int minK, int minLegitimateSpacer, int maxLegitimateSpacer, unordered_map<string,double>& stats, bool strict, 
                     bool preStrict, ofstream& logFile, int interval, int maxK);

bool skipThisLine(const string& read, double iligitimateRatio);

void findSeedPattern(string line, unordered_map<string,vector<int>>& singleLineMapSeedKToIdx, int seedK);

bool notOverlapping(int idxStartPotential, int idxStartCompare, int idxEndPotential, int idxEndCompare, 
                    int minLegitimateSpacer, int maxLegitimateSpacer);
void expandSeedToKmer(const string& line, string smer, vector<int> smerIdxVect , int minK, 
                    set<string>& uniqueKmersInLine, int minLegitimateSpacer, int maxLegitimateSpacer, bool strict, int maxK);

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

unordered_map<int, string> reCannonization(const unordered_map<int, string>& provisionalRepList, const DynamicBins& bins, ofstream& logFile);

void creatingOutputMap(unordered_map<string,data_t>& outputMap, unordered_map<string,data_t>& binsOutputMap, 
                       const unordered_map<int,string>& choosenReps, const unordered_map<int,
                       vector<string>>& reverseBins, const unordered_map<string,data_t>& kmap);

int palindromicScore(string kmer, int alpha = 0);

void validateBins(const unordered_map<int, string>& provisionalRepList, const DynamicBins& bins, 
                  const unordered_map<int,vector<string>>& reverseBins, ofstream& logFile);

// bool isKmer(string maybeKmer);

// step three functions:
int buildSmap(ifstream& inCatalog, unordered_map<string,Kmap_t>& smap, int seedK);

void findKmersInFileWithSmap(MultiFormatFileReader& fileReader, unordered_map<string,data_t>& globalKmerMap, 
                            unordered_map<string,Kmap_t>& smap, int seedK, unordered_map<string,double>& stats, ofstream& logFile, int minLegitimateSpacer, int minK, int interval);
void expandSeedToKmerWithSmap(const string& line, const string& smer, int& idxInLine, unordered_map<string,data_t>& globalKmerMap, unordered_map<string,Kmap_t>& smap, bool& activeLine, unordered_map<string,int>& kmerToIdxInLine);

bool willSelfOverlap(const unordered_map<string,int>& kmerToIdxInLine,int startIdexOfKmerInLine, string kmerInLine);

bool valideHeader(string header);

// Array Dump functions:
void arrayIdentifior(MultiFormatFileReader& fileReader, 
                     unordered_map<string,data_t>& globalArrayMap, 
                     unordered_map<string,Kmap_t>& smap, 
                     int seedK, 
                     unordered_map<string,double>& stats, 
                     ofstream& logFile, 
                     int minLegitimateSpacer, 
                     int maxLegitimateSpacer, 
                     int minK, 
                     int interval);
string expandSeedToKmer(const string& line, const string& smer, int& idxInLine, unordered_map<string,Kmap_t>& smap, bool& activeLine, int& tempStartIdx);



#endif /* functions_hpp */
