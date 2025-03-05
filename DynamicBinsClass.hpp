//
//  DynamicBinsClass.hpp
//  repeatesSearchProjectAltK
//
//  Created by Sarah Katz on 2/3/25.
//

#ifndef DynamicBinsClass_hpp
#define DynamicBinsClass_hpp

#include <stdio.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include "classData.hpp"
#include "globalFunctions.h"

using namespace std;

class DynamicBins {
public:
    DynamicBins();
    void addAutoSingle(const string &Kmer);
    void autoGroup(const string &Kmer1, const string &Kmer2);
    void addToExistingBin(const string &Kmer, int binNumber);
    int getBin(const string &Kmer) const;
    void merge(int binA , int binB);
    set<int> getAssignedBinNumbers() const;
    void print() const;
    unordered_map<string, int> getBinsStructure() const;
    void reverseBin();
    void normlizeReverseBin();
    unordered_map<int,vector<string>> getReBins();
    // void iterReBins(const unordered_map<string,int>& Kmap, function<void(const vector<string>&, int, const unordered_map<string,int>&)> func) const;
    int getLen();
    int getNumBins();

private:
    unordered_map<string, int> bins;
    unordered_map<int, vector<string>> reBins;
    int nextBinNumber;
    set<int> existingBinNums; // record of bin nums in use
};

#endif /* DynamicBinsClass_hpp */
