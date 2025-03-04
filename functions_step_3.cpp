#include "functions.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"


// this function initilizes and builds and Smap that maps from Smers to 
// and the indecies / positions where they appear in in the Kmers
void buildSmap(ifstream& InCatalog, unordered_map<string,Kmap_t>& Smap, int seedK) {
    // create a variable for a temp line and a Kmer to extract the Kmer from the catalog output
    string tempLine;
    string Kmer;
    // dump the header line
    getline(InCatalog, tempLine);
    while (!InCatalog.eof() && InCatalog.good()) {
        // while the file is good and doesnt reach its end we keep pulling lines and extracting the first string (the Kmer)
        getline(InCatalog, tempLine);
        istringstream iss(tempLine);
        if (iss >> Kmer) {
            // we also want to extract the reverse of the Kmer so that we can search for both
            string reverseComp = reverseComplement(Kmer);
            // we iterate over the Kmer or reverse complemet 
            // untill our iterator reaches a length that couldn't produce an Smer (end of line - seed)
            for (int j = 0; j <= (Kmer.size() - seedK); j++){
                // we generate a "forward" Smer of size seedK at every index
                string Smer = Kmer.substr(j,seedK);
                auto& Kmap = Smap[Smer];
                auto& idxVect = Kmap[Kmer];
                idxVect.emplace_back(j); // we store the position of the Smer in the Kmer
                
                // we generate a reverse Smer of size seedK at every index
                string reverseSmer = reverseComp.substr(j,seedK);
                auto& KmapReverse = Smap[reverseSmer];
                auto& idxVectReverse = KmapReverse[reverseComp];
                idxVectReverse.emplace_back(j); // we store the position of the Smer in the Kmer
            }
        }
    }
}

// this function goes line by line and searches for known Smers to find known Kmers
void findKmersInFileWithSmap(MultiFormatFileReader& fileReader, unordered_map<string,data_t>& globalKmerMap, unordered_map<string,Kmap_t>& Smap, int seedK){
    // we initilize a line variable and reassign the next line to it using the filereader class untill the file's end
    string line;
    while (fileReader.getNextLine(line)) {
        // we iterate over the line (as long as not empty), generating an Smer at each index point 
        // untill the point we'd pass the end of the line if we made another Smer
        if (line.size() == 0){
            continue;
        }
        for (int i = 0; i <= (line.size() - seedK); i++){
            string Smer = line.substr(i,seedK);
            // if this Smer appears in our Smap we try expanding it to a Kmer and checking if its a known Kmer
            if (Smap.find(Smer) != Smap.end()){
                expandSeedToKmerWithSmap(line, Smer, i, globalKmerMap, Smap);
                // if this succeeds the Kmer is counted in or added to the global Kmer map
            }
        }
        // We iterate over entire global Kmap to 0 out the count in line since we are moving to the next line
        for (const auto& [Kmer, data] : globalKmerMap){
            auto& finalData = globalKmerMap.at(Kmer);
            finalData.countInLine = 0;
        }
    }
}

// this function checks if the known Smer occurance indeed means a known Kmer occurance and if 
void expandSeedToKmerWithSmap(const string& line, const string& Smer, int& idxInLine, unordered_map<string,data_t>& globalKmerMap, unordered_map<string,Kmap_t>& Smap){
    // we acess the Kmers that the Smer of interest appears in
    auto& Kmap = Smap[Smer];
    int copyIdxInLine = idxInLine;
    // we iterate over all the Kmers that this Smer is associated to
    for (const auto& [Kmer, idxs] : Kmap){
        // we generate an empty string for the Kmer in line (that may or may not be our Kmer of interest)
        string KmerInLine;
        // we record the recorded Kmer's length
        int K = Kmer.size();
        // then we iterate over all the indecies that the Smer of interest appears at in this particular recorded Kmer 
        for (int i = 0; i < idxs.size(); i++){
            // an index where this Smer appears in the Kmer
            int startIdxInKmer = idxs[i];
            // a check to prevent substring out of bounds error
            if (copyIdxInLine >= startIdxInKmer){
                // we extract the substring of the potential Kmer from the line
                KmerInLine = line.substr((copyIdxInLine - startIdxInKmer), K);
            }
            else{
                continue;
            }
            // if the Kmer that we extracted is the same Kmer as the recorded Kmer:
            // 1) we cannonize it so that we dont record both Kmer and reverse complement seperatly
            // 2) we count it in the file and in the line
            // 3) we test if its only been counted once in the line, and if so we also count a line it appeared on
            // 4) last we move the iterator of the line in the external function that iterates over the line to pass the Kmer
            if (Kmer == KmerInLine){
                string canonizedKmer = pickKey(Kmer);
                auto& finalData = globalKmerMap[canonizedKmer];
                finalData.countInFile++;
                finalData.countInLine++;
                if (finalData.countInLine <= 1){
                    finalData.numLines++;
                }
                // to move iterator we:
                // check that its still the original index that the external function gave, if so we update
                if (idxInLine == copyIdxInLine){
                    idxInLine += (Kmer.length() - startIdxInKmer);
                }
                // if it isn't, then we chack if the previous K that we found at this Smer was smaller then this one
                // if so, we update to what the skip over this K would have been
                else if(idxInLine - copyIdxInLine < Kmer.length()){
                    idxInLine = copyIdxInLine + (Kmer.length() - startIdxInKmer);
                }
            }
        }
    }
}