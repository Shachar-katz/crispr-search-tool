#include "DynamicBinsClass.hpp"

using namespace std;

DynamicBins::DynamicBins() : nextBinNumber(1) {}

    // adds a key with an auto-generated bin number
void DynamicBins::addAutoSingle(const string &kmer) {
    if(this->getBin(kmer) != -1){
        return;
    }
    int binNumber = nextBinNumber++;
    bins[kmer] = binNumber;
    existingBinNums.insert(binNumber);
}
    // groups kmers with an auto generated bin number
void DynamicBins::autoGroup(const string &kmer1, const string &kmer2) {
    int binNumber = nextBinNumber++;
    bins[kmer1] = binNumber;
    bins[kmer2] = binNumber;
    existingBinNums.insert(binNumber);
}
    // Manually assigns a bin number to a key.
void DynamicBins::addToExistingBin(const string &kmer, int binNumber) {
    bins[kmer] = binNumber;
}

    // Looks up the bin number for a key or -1 if not found
int DynamicBins::getBin(const string &kmer) const {
    if (bins.find(kmer) != bins.end()){
        return bins.at(kmer);
    }
    else{
        return -1;
    }
}

    // mergers bins to the smaller bin 
void DynamicBins::merge(int binA, int binB) {
    if(binA == binB){
        return;
    } // viable_change
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
    for (const auto& [kmer, binNum] : bins) {
        auto& KVect = reBins[binNum];
        KVect.emplace_back(kmer);
    }
    normlizeReverseBin();
}
void DynamicBins::normlizeReverseBin(){
    for (auto& [binNum, KVect] : reBins) {
        sort(KVect.begin(), KVect.end(),[](const std::string &a, const std::string &b) {
                return (a.length() < b.length()) || (a.length() == b.length() && pickKey(a) < pickKey(b));
            }
        );
    }
}

unordered_map<int,vector<string>> DynamicBins::getReBins(){
    return reBins;
}

    // gets len of bins structure
int DynamicBins::getLen(){
    return bins.size();
}

int DynamicBins::getNumBins(){
    return this->existingBinNums.size();
}
