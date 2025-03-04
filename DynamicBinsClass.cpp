#include "DynamicBinsClass.hpp"
#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;

DynamicBins::DynamicBins() : nextBinNumber(1) {}

    // adds a key with an auto-generated bin number
void DynamicBins::addAutoSingle(const string &Kmer) {
    if(this->getBin(Kmer) != -1){
        return;
    }
    int binNumber = nextBinNumber++;
    bins[Kmer] = binNumber;
    existingBinNums.insert(binNumber);
}
    // groups Kmers with an auto generated bin number
void DynamicBins::autoGroup(const string &Kmer1, const string &Kmer2) {
    int binNumber = nextBinNumber++;
    bins[Kmer1] = binNumber;
    bins[Kmer2] = binNumber;
    existingBinNums.insert(binNumber);
}
    // Manually assigns a bin number to a key.
void DynamicBins::addToExistingBin(const string &Kmer, int binNumber) {
    bins[Kmer] = binNumber;
}

    // Looks up the bin number for a key or -1 if not found
int DynamicBins::getBin(const string &Kmer) const {
    if (bins.find(Kmer) != bins.end()){
        return bins.at(Kmer);
    }
    else{
        return -1;
    }
}

    // mergers bins to the smaller bin 
void DynamicBins::merge(int binA, int binB) {
    int mergedBin = std::min(binA, binB);
    int removedBin = std::max(binA, binB);
    
    for (auto & pair : bins) {
        if (pair.second == removedBin) {
            pair.second = mergedBin;
        }
    }
    existingBinNums.erase(removedBin);
    existingBinNums.insert(mergedBin);
}

    // Returns all currently assigned bin numbers.
set<int> DynamicBins:: getAssignedBinNumbers() const {
        return existingBinNums;
}

    // Prints the mapping from keys to bin numbers for debugging
void DynamicBins::print() const {
        for (const auto &pair : bins)
            cout << "Key: " << pair.first << " -> Bin number: " << pair.second << "\n";
}

void DynamicBins::reverseBin(){
    for (const auto& [Kmer, binNum] : bins) {
        auto& KVect = reBins[binNum];
        KVect.emplace_back(Kmer);
    }
}

unordered_map<int,vector<string>> DynamicBins::getReBins(){
    return reBins;
}

//     // iterates over the bins the reversal - NOT NEEDED 
// void DynamicBins::iterReBins(const unordered_map<string,int>& Kmap, function<void(const vector<string>&, int, const unordered_map<string,int>&)> func) const {
//     for (const auto& [binNum, KVect] : reBins) {
//         try {
//             func(KVect, binNum, Kmap);
//         } catch (out_of_range) {
//             cerr << "here is caught for key:" << Kmer << endl;
//             throw -1;
//         }
        
//     }
// }

    // gets len of bins structure
int DynamicBins::getLen(){
    return bins.size();
}

int DynamicBins::getNumBins(){
    return this->existingBinNums.size();
}
