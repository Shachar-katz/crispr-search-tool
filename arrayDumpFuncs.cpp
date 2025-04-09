#include "functions.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"
#include "arrayClass.h"

// this function goes line by line and searches for arrays by known repeat identification
void arrayIdentifior(MultiFormatFileReader& fileReader, 
                     unordered_map<string,data_t>& globalArrayMap, 
                     unordered_map<string,Kmap_t>& smap, 
                     int seedK, 
                     unordered_map<string,double>& stats, 
                     ofstream& logFile, 
                     int minLegitimateSpacer, 
                     int maxLegitimateSpacer, 
                     int minK, 
                     int interval)
{
    // stats
    int progressCounter = 0;
    int numReadsWithArrays = 0; // ?? potentially add numReads with array from a particular repeat
    // line variable and filereader class used to iterate over lines in the file
    string line;
    while (fileReader.getNextLine(line)) {
        bool activeLine = false; // stats
        LineArrayHandler arrayHandler(line, maxLegitimateSpacer);        
        // untill the point we'd pass the end of the line if we made another smer
        if (line.length() <= (2 * minK + minLegitimateSpacer + 2)){ continue; }
        // iterate over the line in jumps of s, and generate an smer at each index point 
        int i = 0;
        while (i <= (line.size() - seedK)){
            string smer = line.substr(i,seedK);
            // if this smer appears in our smap we try expanding it to a kmer and checking if its a known kmer
            if (smap.find(smer) == smap.end()){ 
                int tempStartIdx;
                string repeat = expandSeedToKmer(line, smer, i, smap, activeLine, tempStartIdx);
                i++; // for now to keep everything stable !!
                if (repeat == "") { continue; }
                arrayHandler.manageState(repeat, tempStartIdx);
            }
            else{
                i += seedK;
            }
        }
        if (arrayHandler.isActive()) { arrayHandler.uploadArray(); }
        //@Here ?? now what do i do with the line array vect? 
        // stats
        if (activeLine == true){ numReadsWithArrays++; }
        // progress bar
        progressCounter++;
        if (progressCounter % interval == 0){
            cout << "Procession line: " << progressCounter << endl;
            logFile << "Procession line: " << progressCounter << endl;
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

// this function checks if the known smer occurance indeed means a known kmer occurance and if 
string expandSeedToKmer(const string& line, const string& smer, int& idxInLine, unordered_map<string,Kmap_t>& smap, bool& activeLine, int& tempStartIdx){
    // we acess the Kmers that the smer of interest appears in
    auto& kmap = smap[smer];
    // we iterate over all the Kmers that this smer is associated to
    for (const auto& [kmer, idxs] : kmap){
        // we generate an empty string for the kmer in line (that may or may not be our kmer of interest)
        string kmerInLine;
        // we record the recorded kmer's length
        int k = kmer.size();
        // iterate over all the indecies that the smer of interest appears at in this particular recorded kmer 
        for (int i = 0; i < idxs.size(); i++){
            // an index where this smer appears in the kmer
            int startIdxInKmer = idxs[i];
            int startIdexOfKmerInLine = idxInLine - startIdxInKmer;
            // a check to prevent substring out of bounds error
            if (idxInLine < startIdxInKmer){ continue; }
            // extract the substring of the potential kmer from the line
            kmerInLine = line.substr(startIdexOfKmerInLine, k);
            if (kmer != kmerInLine){ continue; } // early rejection if kmer not a match

            activeLine = true;
            idxInLine += (kmer.length() - startIdxInKmer);
            tempStartIdx = startIdexOfKmerInLine;
            return kmer;
        }
    }
    return "";
}