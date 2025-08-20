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
                     unordered_map<string,ArrayPositionData>& arrayPositionMap, 
                     unordered_map<string,Kmap_t>& smap,
                     unordered_map<string,string>& kmerToId, 
                     int seedK, 
                     unordered_map<string,double>& stats, 
                     ofstream& logFile,
                     ofstream& readDump,
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
                int spacerLen = 0;
                if (activeLine) { spacerLen = arrayHandler.getSpacerLen(); }
                string repeat = expandSeedToKmer(line, smer, i, spacerLen, smap, activeLine, tempStartIdx, maxMismatches, numMissmatches);
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
        readDump << "> R_" << progressCounter << endl << line << endl;
        vector<Array> lineArrayVect = arrayHandler.getLineArrayVect();
        for (auto& Array : lineArrayVect){
            string arrId = "A_" + to_string(arrayId);
            globalArrayMap[arrId] = Array;

            string readId = "R_" + to_string(progressCounter);
            ArrayPositionData currData;
            currData.readID = readId;
            currData.startPos = Array.getStartPos();
            currData.endPos = Array.getEndPos();
            arrayPositionMap[arrId] = currData;
            
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


// this function checks if the known smer occurance means a known kmer occurance 
string expandSeedToKmer(const string& line, 
                        const string& smer, 
                        int& idxInLine, 
                        const int& spacerLen,
                        unordered_map<string,Kmap_t>& smap, 
                        bool& activeLine, 
                        int& tempStartIdx, 
                        int maxMismatches, 
                        int& numMissmatches){
    // we acess the Kmers that the smer of interest appears in
    auto& kmap = smap[smer];
    // Create a vector of kmers and sort by size (largest first)
    vector<string> sortedKmers;
    sortedKmers.reserve(kmap.size()); // Reserve space for efficiency
    
    for (const auto& [kmer, idxs] : kmap) {
        sortedKmers.push_back(kmer);
    }
    
    // Sort by size (largest first), with tie-breaking by canonical form for consistency
    sort(sortedKmers.begin(), sortedKmers.end(), [](const string& a, const string& b) {
        if (a.size() != b.size()) {
            return a.size() > b.size(); // Sort by size, largest first
        }
        return pickKey(a) < pickKey(b); // Tie-breaker using canonical form
    });
    // we iterate over all the Kmers that this smer is associated to
    for (const string& kmer : sortedKmers) {
        const auto& idxs = kmap[kmer];
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
            if (!isKmerMatch(line, startIdxOfKmerInLine, endIdxOfKmerInLine, smerSet, smer.size(), numMissmatches, maxMismatches)) { break; }

            activeLine = true;
            idxInLine += (kmer.length() - startIdxInKmer + (spacerLen)); // update index
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