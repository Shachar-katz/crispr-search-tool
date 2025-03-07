//
//  functions.hpp
//  repeatesSearchProject
//
//  Created by Sarah Katz on 10/23/24.
//

#ifndef functions_hpp
#define functions_hpp

#include <stdio.h>
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
bool isInputFileValid(ifstream& inFS);
void findKmersInFile(MultiFormatFileReader& fileReader, unordered_map<string,int>& globalKmerMap, int seedK, int minK, int legitimateSpacer, unordered_map<string,double>& stats);
bool skipThisLine(const string& read);
void findSeedPattern(string line, unordered_map<string,vector<int>>& singleLineMapSeedKToIdx, int seedK);
bool notOverlapping(int idxStartPotential, int idxStartCompare, int idxEndPotential, int idxEndCompare, int legitimateSpacer);
void expandSeedToKmer(const string& line, string Smer, vector<int> SmerIdxVect , int minK, set<string>& uniqueKmersInLine, int legitimateSpacer);

//step two functions:
void catalogToSAndKMaps(ifstream& InCatalog, unordered_map<string,vector<string>>& Smap, unordered_map<string,data_t>& Kmap, int seedK);
void findSmersVect(string shortestKmer, vector<string>& smerVect, int seedK);
void makePotentialRelationsSet(const unordered_map<string,vector<string>>& Smap, set<pair<string,string>>& potentialRelationSet);
void chooseShortestK (const pair<string, string>& Kpair, string& shortestK, string& otherK);
void binRelatives(const string& shortestK, const string& otherK, DynamicBins& bins);
void verifyRelation(const unordered_map<string,vector<string>>& Smap, const set<pair<string,string>>& potentialRelationSet, int seedK, DynamicBins& bins, int alpha = 1);
string kmerCompetition(const unordered_map<string,data_t>& Kmap, string currentRep, string auditioningKmer);
string tieBreaker(string currentRep, string auditioningKmer);
void binSingles(const unordered_map<string,vector<string>>& Smap, DynamicBins& bins);
void selectReps(unordered_map<int, string>& provisionalRepList, const unordered_map<int,vector<string>>& reverseBins, const unordered_map<string,data_t>& Kmap);
unordered_map<int, string> reCannonization(const unordered_map<int, string>& provisionalRepList, const DynamicBins& bins);
void creatingOutputMap(unordered_map<string,data_t>& outputMap, unordered_map<string,data_t>& binsOutputMap, const unordered_map<int,string>& choosenReps, const unordered_map<int,vector<string>>& reverseBins, const unordered_map<string,data_t>& Kmap);
int palindromicScore(string Kmer);

// step three functions:
void buildSmap(ifstream& InCatalog, unordered_map<string,Kmap_t>& Smap, int seedK);
void findKmersInFileWithSmap(MultiFormatFileReader& fileReader, unordered_map<string,data_t>& globalKmerMap, unordered_map<string,Kmap_t>& Smap, int seedK, unordered_map<string,double>& stats);
void expandSeedToKmerWithSmap(const string& line, const string& Smer, int& idx , unordered_map<string,data_t>& globalKmerMap, unordered_map<string,Kmap_t>& Smap, bool& activeLine);


#endif /* functions_hpp */
