#include "DynamicBinsClass.hpp"

using namespace std;

DynamicBins::DynamicBins(int maxClusterSize, int seedK) 
    : nextBinNumber(1), 
      maxClusterSize(maxClusterSize), seedK(seedK) {}

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

int DynamicBins::getNextBinNum() const{
    return nextBinNumber;
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

void DynamicBins::checkAndRecluster(int binNumber) {
    // get bin
   vector<string> clusterKmers = getBinContent(binNumber);
   
   if (clusterKmers.size() > maxClusterSize) {
       // split by GC composition
       vector<vector<string>> primarySubclusters = splitClusterByComposition(clusterKmers);
    //    vector<vector<string>> secondarySubClusters;
       
       // split by complexity if still too large
    //    for(auto& subCluster : primarySubclusters) {
    //        if (subCluster.size() > maxClusterSize) {
    //            vector<vector<string>> complexityGroups = splitByComplexity(subCluster);
    //            for (auto& complexityGroup : complexityGroups) {
    //                secondarySubClusters.push_back(complexityGroup);
    //            }
    //        } else {
    //            secondarySubClusters.push_back(subCluster);
    //        }
    //    }
              
    //    for(auto& subCluster : secondarySubClusters){
    //         cout << subCluster.size() << endl;
    //    } // db
       
       reassignCluster(binNumber, primarySubclusters);
   }
}

void DynamicBins::forceReclusterAll() {
    vector<int> binsToCheck;
    for (int binNum : existingBinNums) {
        binsToCheck.push_back(binNum);
    }
    
    for (int binNum : binsToCheck) {
        checkAndRecluster(binNum);
    }
}

vector<string> DynamicBins::getBinContent(int binNumber) {
    vector<string> kmers;
    for (const auto& [kmer, bin] : bins) {
        if (bin == binNumber) {
            kmers.push_back(kmer);
        }
    }
    return kmers;
}

vector<vector<string>> DynamicBins::splitClusterByComposition(const vector<string>& cluster) {    
    // Simple grouping by GC content ranges
    map<int, vector<string>> gcGroups;
    
    for (const auto& kmer : cluster) {
        int gcCount = 0;
        for (char c : kmer) {
            if (c == 'G' || c == 'C') gcCount++;
        }
        int gcPercent = (gcCount * 100 / kmer.length()) / 20 * 20; // Round to nearest 20%
        gcGroups[gcPercent].push_back(kmer);
    }
    
    vector<vector<string>> result;
    for (auto& [gc, kmers] : gcGroups) {
        result.push_back(kmers);
        cout << "  GC group " << gc << "%: " << kmers.size() << " k-mers" << endl;
    }
    
    return result;
}

vector<vector<string>> DynamicBins::splitByComplexity(const vector<string>& cluster) {
    map<int, vector<string>> complexityGroups;
    
    for (const auto& kmer : cluster) {
        double complexity = calculateComplexity(kmer);
        int complexityBin = (int)(complexity * 4);
        complexityGroups[complexityBin].push_back(kmer);
    }
    
    vector<vector<string>> result;
    for (auto& [complexityLevel, kmers] : complexityGroups) {
        if (!kmers.empty()) {
            result.push_back(kmers);
        }
    }
    
    return result;
}

double DynamicBins::calculateComplexity(const string& kmer) {
    unordered_map<char, int> counts;
    for (char c : kmer) {
        counts[c]++;
    }
    
    double entropy = 0.0;
    int total = kmer.length();
    for (const auto& [nucleotide, count] : counts) {
        if (count > 0) {
            double p = (double)count / total;
            entropy -= p * log2(p);
        }
    }
    
    return entropy;
}

void DynamicBins::reassignCluster(int oldBinNum, const vector<vector<string>>& subclusters) {
    vector<string> oversizedBin = getBinContent(oldBinNum);
    // clear out old bin
    for (const auto& kmer : oversizedBin) {
        bins.erase(kmer);
    }
    existingBinNums.erase(oldBinNum);
    // reassign sublustered kmers
    for (const auto& subcluster : subclusters) {
        if (!subcluster.empty()) {
            int newBinNum = nextBinNumber++;
            for (const auto& kmer : subcluster) {
                bins[kmer] = newBinNum;
            }
            existingBinNums.insert(newBinNum);
        }
    }
}
