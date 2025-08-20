#pragma once
#include <iostream>
#include <string>
#include <unordered_set>
using namespace std;

inline void findSmerSet(string kmer, unordered_set<string>& smerSet, int seedK){
    for (int j = 0; j <= (kmer.size() - seedK); j++){
        string smerForSearch = kmer.substr(j,seedK);
        if (smerForSearch.length() < seedK) {break;}
        smerSet.insert(smerForSearch);
    }
}

inline bool isKmerMatch(const string& line,  
                        int start,
                        int end, 
                        const unordered_set<string>& smerSet, 
                        int seedK, 
                        int& missmatches, 
                        int maxMismatches = 0){
    int dissimilarity = 0;
    
    for (int i = start; i <= (end - seedK); i++){
        if (dissimilarity > maxMismatches) { return false; }
        string smerInLine = line.substr(i,seedK);
        if (smerSet.count(smerInLine) == 0) { dissimilarity++; }
    }
    missmatches = dissimilarity;
    return true;
}

inline bool areRepeatsTheSame(const string& repeatA,
                              const string& repeatB,
                              int seedK, 
                              int maxMismatches = 0){
    int dissimilarity = 0;
    string shortest;
    string longer;
    if (repeatA.length() < repeatB.length()) {
        shortest = repeatA; 
        longer = repeatB;
    }
    else { 
        shortest = repeatB;
        longer = repeatA;
    }
    unordered_set<string> smerSetLonger;
    findSmerSet(longer, smerSetLonger, seedK);
    for (int i = 0; i <= (shortest.length() - seedK); i++){
        if (dissimilarity > maxMismatches) { return false; }
        string smerOfShortest = shortest.substr(i,seedK);
        if (smerSetLonger.count(smerOfShortest) == 0) { dissimilarity++; }
    }
    return true;
}