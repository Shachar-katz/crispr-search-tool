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
    unordered_set<string> smerSetB;
    findSmerSet(repeatB, smerSetB, seedK);
    for (int i = 0; i <= (repeatA.length() - seedK); i++){
        if (dissimilarity > maxMismatches) { return false; }
        string smerInA = repeatA.substr(i,seedK);
        if (smerSetB.count(smerInA) == 0) { dissimilarity++; }
    }
    return true;
}