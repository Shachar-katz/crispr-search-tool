#include "DynamicBinsClass.hpp"

using namespace std;

DynamicBins::DynamicBins(int maxClusterSize, int seedK) 
    : nextBinNumber(1), autoReclusterEnabled(true), 
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
    if (autoReclusterEnabled) {
        checkAndRecluster(binNumber);
    }
}
    // Manually assigns a bin number to a key.
void DynamicBins::addToExistingBin(const string &kmer, int binNumber) {
    bins[kmer] = binNumber;
    if (autoReclusterEnabled) {
        checkAndRecluster(binNumber);
    }
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
    if (autoReclusterEnabled) {
        checkAndRecluster(mergedBin);
    }
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

void DynamicBins::enableAutoReclustering(bool enable) {
    this->autoReclusterEnabled = enable;
}

void DynamicBins::checkAndRecluster(int binNumber) {
    vector<string> clusterKmers = getCurrentClusterKmers(binNumber);
    
    if (clusterKmers.size() > maxClusterSize) {
        cout << "Auto-reclustering bin " << binNumber << " with " << clusterKmers.size() << " k-mers" << endl;
        
        vector<vector<string>> subclusters = splitClusterByComposition(clusterKmers);
        
        if (subclusters.size() > 1) {
            cout << "Split by composition into " << subclusters.size() << " subclusters" << endl;
            reassignCluster(binNumber, subclusters);
            return;
        }
        
        subclusters = splitClusterByConnectivity(clusterKmers);
        if (subclusters.size() > 1) {
            cout << "Split by connectivity into " << subclusters.size() << " subclusters" << endl;
            reassignCluster(binNumber, subclusters);
        }
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

vector<string> DynamicBins::getCurrentClusterKmers(int binNumber) {
    vector<string> kmers;
    for (const auto& [kmer, bin] : bins) {
        if (bin == binNumber) {
            kmers.push_back(kmer);
        }
    }
    return kmers;
}

// vector<vector<string>> DynamicBins::splitClusterByComposition(const vector<string>& cluster) {
//     vector<vector<string>> groups;
//     vector<bool> assigned(cluster.size(), false);
    
//     for (int i = 0; i < cluster.size(); i++) {
//         if (assigned[i]) continue;
        
//         vector<string> group;
//         group.push_back(cluster[i]);
//         assigned[i] = true;
        
//         for (int j = i + 1; j < cluster.size(); j++) {
//             if (assigned[j]) continue;
            
//             if (areCompositionallySimilar(cluster[i], cluster[j])) {
//                 group.push_back(cluster[j]);
//                 assigned[j] = true;
//             }
//         }
        
//         groups.push_back(group);
//     }
    
//     return groups;
// }

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

double DynamicBins::calculateGCContent(const string& kmer) {
    int gcCount = 0;
    for (char c : kmer) {
        if (c == 'G' || c == 'C') gcCount++;
    }
    return (double)gcCount / kmer.length();
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

bool DynamicBins::areCompositionallySimilar(const string& kmer1, const string& kmer2) {
    double gc1 = calculateGCContent(kmer1);
    double gc2 = calculateGCContent(kmer2);
    double complexity1 = calculateComplexity(kmer1);
    double complexity2 = calculateComplexity(kmer2);
    
    const double GC_THRESHOLD = 0.20;
    const double COMPLEXITY_THRESHOLD = 0.5;
    
    bool gcSimilar = abs(gc1 - gc2) <= GC_THRESHOLD;
    bool complexitySimilar = abs(complexity1 - complexity2) <= COMPLEXITY_THRESHOLD;
    
    return gcSimilar && complexitySimilar;
}

vector<vector<string>> DynamicBins::splitClusterByConnectivity(const vector<string>& cluster) {
    unordered_map<string, vector<string>> adjacency;
    
    for (int i = 0; i < cluster.size(); i++) {
        for (int j = i + 1; j < cluster.size(); j++) {
            if (areDirectlySimilar(cluster[i], cluster[j])) {
                adjacency[cluster[i]].push_back(cluster[j]);
                adjacency[cluster[j]].push_back(cluster[i]);
            }
        }
    }
    
    unordered_set<string> visited;
    vector<vector<string>> subclusters;
    
    for (const auto& kmer : cluster) {
        if (visited.count(kmer) == 0) {
            vector<string> component;
            dfsComponent(kmer, adjacency, visited, component);
            subclusters.push_back(component);
        }
    }
    
    return subclusters;
}

bool DynamicBins::areDirectlySimilar(const string& kmer1, const string& kmer2) {
    unordered_set<string> smers1, smers2;
    findSmerSet(kmer1, smers1);
    findSmerSet(kmer2, smers2);
    
    int sharedCount = 0;
    for (const auto& smer : smers1) {
        if (smers2.count(smer)) {
            sharedCount++;
        }
    }
    
    int minSmers = min(smers1.size(), smers2.size());
    int requiredShared = max(1, minSmers / 3);
    
    return sharedCount >= requiredShared;
}

void DynamicBins::findSmerSet(const string& kmer, unordered_set<string>& smerSet) {
    for (int i = 0; i <= (int)kmer.size() - seedK; i++) {
        string smer = kmer.substr(i, seedK);
        smerSet.insert(smer);
    }
}

void DynamicBins::dfsComponent(const string& kmer, 
                              const unordered_map<string, vector<string>>& adjacency,
                              unordered_set<string>& visited, 
                              vector<string>& component) {
    visited.insert(kmer);
    component.push_back(kmer);
    
    auto it = adjacency.find(kmer);
    if (it != adjacency.end()) {
        for (const auto& neighbor : it->second) {
            if (visited.count(neighbor) == 0) {
                dfsComponent(neighbor, adjacency, visited, component);
            }
        }
    }
}

void DynamicBins::reassignCluster(int oldBinNum, const vector<vector<string>>& subclusters) {
    vector<string> oldKmers = getCurrentClusterKmers(oldBinNum);
    for (const auto& kmer : oldKmers) {
        bins.erase(kmer);
    }
    existingBinNums.erase(oldBinNum);
    
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
