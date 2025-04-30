#include "functions.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"
#include "arrayClass.h"

int buildKmap(ifstream& inCatalog, unordered_map<string,string>& kmerToId, int minPalindromic){
    string tempLine;
    // dump the header line as long as valid
    getline(inCatalog, tempLine);
    if (!valideHeader(tempLine)){ return 1; }
    while (!inCatalog.eof() && inCatalog.good()) {
        // iterate over catalog to extract the data
        getline(inCatalog, tempLine);
        istringstream iss(tempLine);
        string kmer; // we keep
        int number_of_lines;
        string binId; // we keep
        int palindromic;
        int lengthK;
        // Attempt to parse
        if (!(iss >> kmer >> number_of_lines >> binId >> palindromic >> lengthK) && !tempLine.empty()) {
            cerr << "Could not parse line: " << tempLine << endl;
            continue;
        }
        if (palindromic < minPalindromic){ continue; }
        // we also want to extract the reverse of the kmer so that we can search for both
        string reverseComp = reverseComplement(kmer);
        kmerToId[kmer]    = binId + "_a";
        kmerToId[reverseComp] = binId + "_b";
    }
    return 0;
}


// this function goes line by line and searches for arrays by known repeat identification
void arrayIdentifior(MultiFormatFileReader& fileReader, 
                     unordered_map<string,Array>& globalArrayMap, 
                     unordered_map<string,Kmap_t>& smap,
                     unordered_map<string,string>& kmerToId, 
                     int seedK, 
                     unordered_map<string,double>& stats, 
                     ofstream& logFile, 
                     int minLegitimateSpacer, 
                     int maxLegitimateSpacer, 
                     int minK, 
                     int interval,
                     int maxMismatches)
{
    // stats
    int progressCounter = 0;
    int numReadsWithArrays = 0; // ?? potentially add numReads with array from a particular repeat
    int arrayId = 1;
    string line;
    while (fileReader.getNextLine(line)) {
        bool activeLine = false; // stats
        // in line array handler object
        LineArrayHandler arrayHandler(line, maxLegitimateSpacer, minLegitimateSpacer);        
        // if the line is too short we skip it
        if (line.length() <= (2 * minK + minLegitimateSpacer + 2)){ continue; }
        // iterate over the line in jumps of s, and generate an smer at each index point 
        int i = 0;
        while (i <= (line.size() - seedK)){
            string smer = line.substr(i,seedK);
            // if this smer appears in our smap we try expanding it to a kmer and checking if its a known kmer
            if (smap.find(smer) != smap.end()){ 
                int tempStartIdx;
                int numMissmatches;
                string repeat = expandSeedToKmer(line, smer, i, smap, activeLine, tempStartIdx, maxMismatches, numMissmatches);
                i++; // for now to keep everything stable !! (expand seed to Kmer also skips k nucleotides)
                if (repeat == "") { continue; }
                // if a known Kmer was found:
                // a) search for its ID
                // b) feed it to the array handler to either:
                    // i) expand an existing array 
                    // ii) create a new array
                    // iii) close the existing array and reopen a new array
                string id = kmerToId.at(repeat);
                arrayHandler.manageState(repeat, tempStartIdx, id, numMissmatches);
            }
            else{
                i += seedK;
            }
        }
        // if at the end of the line an array is left open close it
        if (arrayHandler.isActive()) { arrayHandler.uploadArray(); }
        // stats
        if (activeLine == true){ numReadsWithArrays++; }
        // progress bar
        progressCounter++;
        if (progressCounter % interval == 0){
            cout << "Procession line: " << progressCounter << endl;
            logFile << "Procession line: " << progressCounter << endl;
        }
        // upload line to global array map
        if (arrayHandler.noArrays()) { continue; }
        vector<Array> lineArrayVect = arrayHandler.getLineArrayVect();
        for (auto& Array : lineArrayVect){
            string arrId = "A_" + to_string(arrayId);
            globalArrayMap[arrId] = Array;
            arrayId++;
        }
    }
    // check for 0 division (process terminated before it started)
    if (progressCounter == 0){
        cerr << "no lines were processed" << endl;
        return;
    }
    // calculate stats
    double precentReadsWithRepeat = (static_cast<double>(numReadsWithArrays) / static_cast<double>(progressCounter)) * 100;
    stats["number_of_reads_in_file: "] += progressCounter;
    stats["number_of_reads_in_file_with_array: "] += numReadsWithArrays;
    stats["precent_reads_in_file_with_array: "] += precentReadsWithRepeat;
}

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
                        int maxMismatches){
    int dissimilarity = 0;
    
    for (int i = start; i <= (end - seedK); i++){
        if (dissimilarity > maxMismatches) { return false; }
        string smerInLine = line.substr(i,seedK);
        if (smerSet.count(smerInLine) == 0) { dissimilarity++; }
    }
    missmatches = dissimilarity;
    return true;
}

// this function checks if the known smer occurance means a known kmer occurance 
string expandSeedToKmer(const string& line, 
                        const string& smer, 
                        int& idxInLine, 
                        unordered_map<string,Kmap_t>& smap, 
                        bool& activeLine, 
                        int& tempStartIdx, 
                        int maxMismatches, 
                        int& numMissmatches){
    // we acess the Kmers that the smer of interest appears in
    auto& kmap = smap[smer];
    // we iterate over all the Kmers that this smer is associated to
    for (const auto& [kmer, idxs] : kmap){
        // we record the recorded kmer's length
        int k = kmer.size();
        unordered_set<string> smerSet;
        findSmerSet(kmer, smerSet, smer.size());
        // iterate over all the indecies that the smer of interest appears at in this particular recorded kmer 
        for (int i = 0; i < idxs.size(); i++){
            // an index where this smer appears in the kmer
            int startIdxInKmer = idxs[i];
            int startIdxOfKmerInLine = idxInLine - startIdxInKmer;
            int endIdxOfKmerInLine = startIdxOfKmerInLine + k;
            // a check to prevent out of bounds repeat
            if (idxInLine < startIdxInKmer){ continue; }
            if (endIdxOfKmerInLine > line.size()) { continue; }
            // kmerInLine = line.substr(startIdexOfKmerInLine, k);
            // if (kmer != kmerInLine){ continue; } // early rejection if kmer not a match
            maxMismatches = smer.length();
            if (!isKmerMatch(line, startIdxOfKmerInLine, endIdxOfKmerInLine, smerSet, smer.size(), numMissmatches, maxMismatches)) { break; }

            activeLine = true;
            idxInLine += (kmer.length() - startIdxInKmer); // update index
            tempStartIdx = startIdxOfKmerInLine; // update start idx for array
            return kmer; //@here - should i just have it return all of this and log the missmatches and not change the array handler and just have it cut the sections from the line the same way?
        }
    }
    return "";
}

void constructRepeatMap(const unordered_map<string,Array>& globalArrayMap, 
                        unordered_map<string,RepeatData>& repeatMap){
    for (auto& [arrId, array] : globalArrayMap){
        vector<int> tempErrorVect = array.getErrorVect();
        for (int i = 0; i < tempErrorVect.size(); i++){
            RepeatData data;
            data.arrayID = arrId + "_" + to_string(i);
            data.repeatID = array.getRepeatId();
            data.repeatIdx = i;
            data.numMissmatches = tempErrorVect[i];
            data.repeatLen = array.getRepeat().size();
            repeatMap[data.arrayID] = data;
        }
    }
}