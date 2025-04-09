#include "functions.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"
#include "arrayClass.h"
//@here
// this function goes line by line and searches for known Smers to find known Kmers, to find arrays
void findArraysInFileWithSmap(MultiFormatFileReader& fileReader, 
                             unordered_map<string,data_t>& globalArrayMap, 
                             unordered_map<string,Kmap_t>& smap, 
                             int seedK, 
                             unordered_map<string,double>& stats, 
                             ofstream& logFile, 
                             int minLegitimateSpacer, 
                             int minK, 
                             int interval){
    // stats
    int progressCounter = 0;
    int numReadsWithArrays = 0;
    // we initilize a line variable and reassign the next line to it using the filereader class untill the file's end
    string line;
    while (fileReader.getNextLine(line)) {
        // initilize a var to see if this line has a repeat in it (for stats)
        bool activeLine = false;
        unordered_map<string,int> kmerToIdxInLine;
        // we iterate over the line (as long as not empty), generating an smer at each index point 
        // untill the point we'd pass the end of the line if we made another smer
        if (line.length() <= (2 * minK + minLegitimateSpacer + 2)){
            continue;
        }
        int i = 0;
        while (i <= (line.size() - seedK)){
            // cout << "start iter: " << i << endl; // debugg
            string smer = line.substr(i,seedK);
            // if this smer appears in our smap we try expanding it to a kmer and checking if its a known kmer
            if (smap.find(smer) != smap.end()){
                string repeat = arrayExpandSeedToKmer(line, smer, i, smap, activeLine, kmerToIdxInLine);
                i++; // for now to keep everything stable !!
                if (repeat == "") { continue; }

            }
            else{
                i += minK - seedK;
            }
            // cout << "end iter: " << i << endl; // debugg
        }

        // stats
        if (activeLine == true){
            numReadsWithArrays++;
        }
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
    // calculate stats:
    double precentReadsWithRepeat = (static_cast<double>(numReadsWithArrays) / static_cast<double>(progressCounter)) * 100;
    stats["number_of_reads_in_file: "] += progressCounter;
    stats["number_of_reads_in_file_with_repeat: "] += numReadsWithArrays;
    stats["precent_reads_in_file_with_repeat: "] += precentReadsWithRepeat;
}

// this function checks if the known smer occurance indeed means a known kmer occurance and if 
string arrayExpandSeedToKmer(const string& line, const string& smer, int& idxInLine, unordered_map<string,Kmap_t>& smap, bool& activeLine, unordered_map<string,int>& kmerToIdxInLine){
    // we acess the Kmers that the smer of interest appears in
    auto& kmap = smap[smer];
    int copyIdxInLine = idxInLine;
    // we iterate over all the Kmers that this smer is associated to
    for (const auto& [kmer, idxs] : kmap){
        // we generate an empty string for the kmer in line (that may or may not be our kmer of interest)
        string kmerInLine;
        // we record the recorded kmer's length
        int k = kmer.size();
        // then we iterate over all the indecies that the smer of interest appears at in this particular recorded kmer 
        for (int i = 0; i < idxs.size(); i++){
            // an index where this smer appears in the kmer
            int startIdxInKmer = idxs[i];
            int startIdexOfKmerInLine = copyIdxInLine - startIdxInKmer;
            // a check to prevent substring out of bounds error
            if (copyIdxInLine < startIdxInKmer){ continue; }
            // we extract the substring of the potential kmer from the line
            kmerInLine = line.substr(startIdexOfKmerInLine, k);
            if (kmer != kmerInLine){ continue; }
            // if the kmer that we extracted is the same kmer as the recorded kmer:
            activeLine = true;
            // we move the iterator of the line in the external function that iterates over the line to pass the kmer
            // to move iterator we:
            // check that its still the original index that the external function gave, if so we update
            if (idxInLine == copyIdxInLine){
                idxInLine += (kmer.length() - startIdxInKmer);
                kmerToIdxInLine[kmer] = startIdexOfKmerInLine + k - 1;
            }
            // if it isn't, then we chack if the previous k that we found at this smer was smaller then this one
            // if so, we update to what the skip over this k would have been
            else if(idxInLine - copyIdxInLine < kmer.length()){
                idxInLine = copyIdxInLine + (kmer.length() - startIdxInKmer);
                kmerToIdxInLine[kmer] = startIdexOfKmerInLine + k - 1;
            }
            return kmer;
        }
    }
    return "";
}

string expandKmerToArray(const string& line, const string& smer, int& idxInLine, unordered_map<string,Kmap_t>& smap, bool& activeLine, unordered_map<string,int>& kmerToIdxInLine){
    // we acess the Kmers that the smer of interest appears in
    auto& kmap = smap[smer];
    int copyIdxInLine = idxInLine;
    // we iterate over all the Kmers that this smer is associated to
    for (const auto& [kmer, idxs] : kmap){
        // we generate an empty string for the kmer in line (that may or may not be our kmer of interest)
        string kmerInLine;
        // we record the recorded kmer's length
        int k = kmer.size();
        // then we iterate over all the indecies that the smer of interest appears at in this particular recorded kmer 
        for (int i = 0; i < idxs.size(); i++){
            // an index where this smer appears in the kmer
            int startIdxInKmer = idxs[i];
            int startIdexOfKmerInLine = copyIdxInLine - startIdxInKmer;
            // a check to prevent substring out of bounds error
            if (copyIdxInLine < startIdxInKmer){ continue; }
            // we extract the substring of the potential kmer from the line
            kmerInLine = line.substr(startIdexOfKmerInLine, k);
            if (kmer != kmerInLine){ continue; }
            // if the kmer that we extracted is the same kmer as the recorded kmer:
            activeLine = true;
            // we move the iterator of the line in the external function that iterates over the line to pass the kmer
            // to move iterator we:
            // check that its still the original index that the external function gave, if so we update
            if (idxInLine == copyIdxInLine){
                idxInLine += (kmer.length() - startIdxInKmer);
                kmerToIdxInLine[kmer] = startIdexOfKmerInLine + k - 1;
            }
            // if it isn't, then we chack if the previous k that we found at this smer was smaller then this one
            // if so, we update to what the skip over this k would have been
            else if(idxInLine - copyIdxInLine < kmer.length()){
                idxInLine = copyIdxInLine + (kmer.length() - startIdxInKmer);
                kmerToIdxInLine[kmer] = startIdexOfKmerInLine + k - 1;
            }
            return kmer;
        }
    }
    return "";
}