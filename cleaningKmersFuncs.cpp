
#include "functions.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"

// this function populates 2 maps: 
// 1) smap: maps from smer to a vector of Kmers
// 2) kmap: maps from kmer to abundance data that is recived in the file
void catalogToSAndKMaps(ifstream& InCatalog, 
                        unordered_map<string,vector<string>>& smap, 
                        unordered_map<string,data_t>& kmap, 
                        int seedK, 
                        ofstream& logFile) 
{
    logFile << "initializing converting catalog to an smap and kmap" << endl; 
    string line;
    string kmer;
    int abundance;
     // we dump the header
    getline(InCatalog, line);
    // iterate over all lines in the catalog
    while (!InCatalog.eof() && InCatalog.good()) {
        getline(InCatalog, line);
        // extract the line in a stream buffer
        istringstream iss(line);
        // extract the data (kmer string and its abundance)
        if (iss >> kmer >> abundance) {
            if (skipThisLine(kmer, 0.9)){
                continue;
            }
            // we generate a reverse complement for the kmer
            string reverseComp = reverseComplement(kmer);
            // we place the kmer and reverse complement in a map mapping from kmer to abundance, palindrom score, and length
            kmap[kmer].numLines = abundance;
            kmap[kmer].palindromicScore = palindromicScore(kmer);
            kmap[kmer].kLen = kmer.length();
            kmap[reverseComp].numLines = abundance;
            kmap[reverseComp].palindromicScore = palindromicScore(reverseComp);
            kmap[reverseComp].kLen = reverseComp.length();
            // for the length of the kmer - the length of an smer we extract Smers of the kmer and its reverse complement and place them in our smap (maps from smer to vector of Kmers)
            for (int j = 0; j <= (kmer.size() - seedK) ; j++){
                string smer = kmer.substr(j,seedK);
                string reverseSmer = reverseComp.substr(j,seedK);
                if (smer.size() != seedK){ break; }
                smap[smer].emplace_back(kmer);
                smap[reverseSmer].emplace_back(reverseComp);

            }
        }
    }
    logFile << "after initilizing kmap length: " << kmap.size() << "and smap length: " << smap.size() << endl; //debugging
}

// this function recives a kmer input and generates a vector of Smers that comprize it
void findSmersVect(string shortestKmer, vector<string>& smerVect, int seedK){
    for (int j = 0; j <= (shortestKmer.size() - seedK); j++){
        string smerForSearch = shortestKmer.substr(j,seedK);
        // if (SmerForSearch.length() < seedK){ break; }
        smerVect.emplace_back(smerForSearch);
    }
}

// this function populates the potential relations set with pairs of Kmers that are potentially related as seedK of their nucleotides are identical (they are grouped under the same smer)
void makePotentialRelationsSet(const unordered_map<string,vector<string>>& smap, set<pair<string,string>>& potentialRelationSet, ofstream& logFile){
    logFile << "initializing making the potential relations set" << endl; 
    // for every smer and vector of Kmers in the smap, if the vector contains more then 1 kmer:
    for (const auto& [smer, kVect] : smap){
        if (kVect.size() <= 1){ continue; }
        vector<string> copyKVect = kVect;
        // we sort the Kmers from small to large
        sort(copyKVect.begin(), copyKVect.end(), [](const string &a, const string &b) {
            return (a.length() == b.length()) ? pickKey(a) < pickKey(b) : a.length() < b.length();
        });

        // for the sorted vector of Kmers we put every kmer in a pair with all the others to signify they may be subsets of one another
        for (int i = 0; i < copyKVect.size() - 1; i++){
            for (int j = i + 1; j < copyKVect.size(); j++) {
                potentialRelationSet.emplace(copyKVect.at(i), copyKVect.at(j));
            }
        }
    }
    logFile << "after initilizing potential relations set's length: " << potentialRelationSet.size() << endl; //debugging
}

// this is a simple function that gets a pair of strings and 2 empty variables and assigns them them to the variables by length (shortest, and the other one) or if equal chooses order based on cannonic form 
void chooseShortestK (const pair<string, string>& kPair, string& shortestK, string& otherK){
    if (kPair.first.length() < kPair.second.length()) {
        shortestK = kPair.first;
        otherK = kPair.second;
    } else if (kPair.first.length() > kPair.second.length()){
        shortestK = kPair.second;
        otherK = kPair.first;
    } else if (pickKey(kPair.first) < pickKey(kPair.second)){
        shortestK = kPair.first;
        otherK = kPair.second;
    } else{
        shortestK = kPair.second;
        otherK = kPair.first;
    }
}

// this function takes in 2 Kmers that are verified to be related and a bins structure and manages the binning of the Kmers together
void binRelatives(const string& shortestK, const string& otherK, DynamicBins& bins){
    // if the neither Kmers have a bin, we bin them together
    if (bins.getBin(shortestK) == -1 && bins.getBin(otherK) == -1){
        bins.autoGroup(shortestK, otherK);
    }
    // if the shortest k has a bin already, we add the other k to that bin
    else if (bins.getBin(shortestK) != -1 && bins.getBin(otherK) == -1){
        int existingBin = bins.getBin(shortestK);
        bins.addToExistingBin(otherK, existingBin);
    }
    // if the other k has a bin already, we add the shortest k to that bin
    else if (bins.getBin(shortestK) == -1 && bins.getBin(otherK) != -1){
        int existingBin = bins.getBin(otherK);
        bins.addToExistingBin(shortestK, existingBin);
    }
    // if both have bins, we merge the bins
    else if (bins.getBin(shortestK) != -1 && bins.getBin(otherK) != -1){
        int shortestKBin = bins.getBin(shortestK);
        int otherKBin = bins.getBin(otherK);
        bins.merge(shortestKBin, otherKBin);
    }
}

// This function iterates over the pairs of potentially related Kmers and verifies using the smap and a similarity score determined by the size of the Smers and an "alotted error range" variable to determine if they are related 
// If related they are binned together, otherwise they are not.
void verifyRelation(const unordered_map<string,vector<string>>& smap, const set<pair<string,string>>& potentialRelationSet, int seedK, DynamicBins& bins, ofstream& logFile, int alpha){
    logFile << "verifying relations" << endl;
    // for every potentially related kmer pair:
    for (const auto& kPair : potentialRelationSet){
        // we assign them based on length to shorter and longer
        string shortestK;
        string otherK;
        chooseShortestK(kPair, shortestK, otherK);
        // we calculate the required number of Smers that they must NOT appear in together in the vector of Smers of the shorter kmer in order to determine that they are NOT related
        int reqDisimilarity = alpha * seedK; // length diff between 2 k's cant be longer then alpha + 1. ??
        // we then generate a vector of all the Smers that make up the shorter k
        vector<string> smerVect;
        findSmersVect(shortestK, smerVect, seedK);
        int disimilarityScore = 0;
        bool isRelated = true;
        // then we run over this generated smer Vector:
        for(int i = 0; i < smerVect.size(); i++){
            string currentSmer = smerVect[i];
            if (smap.count(currentSmer) == 0){
                logFile << "failed to find " << currentSmer << "in smap" << endl;
                cerr << "failed to find " << currentSmer << "in smap" << endl;
                continue;
            }
            // we then produce for every smer in the smer Vector of our shortest kmer, the vector of Kmers that it appears in on the smap
            const auto& currKvect = smap.at(currentSmer);
            // if we don't find the longer k in this smer we increment the disimilarity score for this pair
            if (find(currKvect.begin(), currKvect.end(), otherK) == currKvect.end()){
                disimilarityScore++;
            }
            // if at any point the disimilarity score passes the required disimilarity score we determine that this pair is NOT related and exit the function
            if (disimilarityScore > reqDisimilarity){
                isRelated = false;
                break;
            }
        }
        // if after the check they are still related, we bin them together using the binning function, as well as their RC's, otherwise we bin them seperately.
        if (isRelated){
            binRelatives(shortestK, otherK, bins);
        }
        else{ 
            bins.addAutoSingle(shortestK);
            bins.addAutoSingle(otherK);
        }
    }
    logFile << "grouping completed" << endl; // debugging
}

// This function iterates over the smap and bins any Kmers that were never put into a pair and thus are not in a bin
void binSingles(const unordered_map<string,vector<string>>& smap, DynamicBins& bins, ofstream& logFile){
    logFile << "begining binning singles" << endl; // debugging
    for (const auto& [smer, kVect] : smap){
        if (kVect.size() == 1){
            string kmer = kVect[0];
            bins.addAutoSingle(kmer);
        }
    }
    logFile << "binning singles" << endl; // debugging
}

// this function takes in 2 Kmers and the kmap and returns the one that is of higher abundance in the kmap
string kmerCompetition(const unordered_map<string,data_t>& kmap, string currentRep, string auditioningKmer, ofstream& logFile){
    // we search for the abundance of both the current leader kmer and the auditioning leader kmer and return the one of higher abundance
    if (kmap.count(currentRep) == 0 || kmap.count(auditioningKmer) == 0){ return ""; }
    int abundanceCurrent = kmap.at(currentRep).numLines;
    int abundanceAudition = kmap.at(auditioningKmer).numLines;

    if (abundanceCurrent < abundanceAudition) {
        return auditioningKmer;
    }
    else if(abundanceCurrent == abundanceAudition){
        return tieBreaker(currentRep,auditioningKmer);
    }
    return currentRep;
}

// this function cannonicly picks a kmer rep if their abundance is equal
string tieBreaker(string currentRep, string auditioningKmer){
    string cannonizedCurrRep = pickKey(currentRep);
    string cannonizedAuditioning = pickKey(auditioningKmer);

    if (cannonizedCurrRep > cannonizedAuditioning){
        return currentRep;
    }
    return auditioningKmer;
}

// this function iterates over the reverse bins and assigns exactly one representative for each bin (the one with the highest abundance).
// the end result is a list of representatives and bin numbers in the repList structure
void selectReps(unordered_map<int, string>& provisionalRepList, const unordered_map<int, vector<string>>& reverseBins, const unordered_map<string,data_t>& kmap, ofstream& logFile){
    logFile << "begining selecting reps" << endl;
    // we iterate over the reverse bins structure that points from bin number to a vector of Kmers in that bin.
    for (const auto& [binNum, kVect] : reverseBins){
        vector <string> sortedKVect = kVect;
        // we sort the Kmers in the vector from small to large
        sort(sortedKVect.begin(), sortedKVect.end(),[](const std::string &a, const std::string &b) {
                return (a.length() < b.length()) || (a.length() == b.length() && pickKey(a) < pickKey(b));
            }
        );
        // for every kmer:
        for (const auto& auditioningKmer : sortedKVect){
            // if that bin doesn't have a representative yet (aka index 0 of our current vector) we assign the kmer as the representative of this bin and continue
            if (provisionalRepList.count(binNum) == 0){
                provisionalRepList[binNum] = auditioningKmer;
                continue;
            }
            // if that bin DOES have a representative and weve reached this part of the code:
            // 1) we assign the representative to a variable and check which is more abundent using kmer competition
            // 2) if the auditioning kmer is more abundent we assign it as the new representative, else we leave the bin as is (with the existing rep)
            string currentRep = provisionalRepList.at(binNum);
            string winner = kmerCompetition(kmap, currentRep, auditioningKmer,logFile);
            if (winner == auditioningKmer){
                provisionalRepList[binNum] = auditioningKmer;
            }
            else if (winner == ""){
                // Check which key is missing
                if (kmap.find(currentRep) == kmap.end()) {
                    logFile << "ERROR: Key not found in kmap: " << currentRep << endl;
                    cerr << "ERROR: Key not found in kmap: " << currentRep << endl;
                }
                else if (kmap.find(auditioningKmer) == kmap.end()) {
                    logFile << "ERROR: Key not found in kmap: " << auditioningKmer << endl;
                    cerr << "ERROR: Key not found in kmap: " << auditioningKmer << endl;
                }
            }
        }
    }
    logFile << "selected final reps" << endl; // debugging
}

// a function to cannonize reps
unordered_map<int, string> reCannonization(const unordered_map<int, string>& provisionalRepList, const DynamicBins& bins, ofstream& logFile) {
    
    unordered_map<string, int> reverseRepList;
    unordered_map<int, string> finalRepList;
    // 1) flip reps list ( kmer -> binNum )
    for(auto& [binNum, repKmer] : provisionalRepList){
        reverseRepList[repKmer] = binNum;
    }
    // 2) iterate over the reverse list, for every kmer create a temp var cannonized = pickkey(kmer)
    for(auto& [repKmer, binNum] : reverseRepList){
        string cannonized = pickKey(repKmer);
        // 3) check if kmer != cannonized && cannonized is another key in reverse rep list:
        if (repKmer != cannonized && (reverseRepList.count(cannonized) != 0)){
            // 4) if so find canonized's binNum, and update that in the new reps list 
            int bin = reverseRepList.at(cannonized);
            finalRepList[bin] = cannonized;
        }
        // 3) otherwise check if the rep is its own cannonized form 
        else if(repKmer == cannonized){
            // 4) if so add it to final reps list
            finalRepList[binNum] = cannonized;
        }
        // 3) otherwise if RC and kmer are in the same bin
        else if(binNum == bins.getBin(cannonized)){
            // 4) if so add it to final reps list
            finalRepList[binNum] = cannonized;
        }
        else{
            logFile << "failed to cannonize the kmer: " << repKmer << endl;
            cerr << "failed to cannonize the kmer: " << repKmer << endl;
        }
    }
    return finalRepList;
}

// this function iterates over all of the structures needed to collect all the data and format it in an output map
void creatingOutputMap(unordered_map<string,data_t>& outputMap, unordered_map<string,data_t>& binsOutputMap, const unordered_map<int,string>& choosenReps, const unordered_map<int, vector<string>>& reverseBins, const unordered_map<string,data_t>& kmap){
    int newBinNum = 1;
    for(const auto& [binNum, repKmer] : choosenReps){
        // filling up normal output
        auto& data = outputMap[repKmer];
        data.binNum = "K_" + to_string(newBinNum);
        data.numLines = kmap.at(repKmer).numLines; 
        data.palindromicScore = kmap.at(repKmer).palindromicScore;
        data.kLen = kmap.at(repKmer).kLen;
        // filling up subKmers output
        auto& subKmersVect = reverseBins.at(binNum); 
        for (string subKmer : subKmersVect){
            auto& dataBins = binsOutputMap[subKmer];
            dataBins.binNum = "K_" + to_string(newBinNum);
            dataBins.palindromicScore = kmap.at(subKmer).palindromicScore;
            dataBins.kLen = kmap.at(subKmer).kLen;
        }
        // updating bin nums
        newBinNum++;
    }
}

// function that calculates palindromic score for repeat (num nucleotides thar are palindromic)
int palindromicScore(string kmer, int alpha){
    int palindromScore = 0;
    int errorCount = 0;
    int iForward = 0;
    int iReverse = kmer.length() - 1;
    while(iForward < iReverse){
        char nucleotide = kmer[iForward];
        char nucleotideRe = kmer[iReverse];
        
        // Directly check if the nucleotides are complementary
        bool isComplementary = false;
        switch(nucleotide) {
            case 'A': isComplementary = (nucleotideRe == 'T'); break;
            case 'G': isComplementary = (nucleotideRe == 'C'); break;
            case 'C': isComplementary = (nucleotideRe == 'G'); break;
            case 'T': isComplementary = (nucleotideRe == 'A'); break;
            default: isComplementary = false;
        }
        
        if (isComplementary) {
            palindromScore++;
        } else if (errorCount < alpha) {
            errorCount++;
        }
        else{
            return palindromScore * 2;
        }
        
        iReverse--;
        iForward++;
    }
    return palindromScore * 2;
}

// function that iteraes over the bins and makes sure they are mirrored
void validateBins(const unordered_map<int, string>& provisionalRepList, const DynamicBins& bins, const unordered_map<int,vector<string>>& reverseBins, ofstream& logFile){
    logFile << "validating bins." << endl;
    bool valid = true;
    for (const auto& [binNum, repK] : provisionalRepList){
        string RC = reverseComplement(repK);
        int binRC = bins.getBin(RC);
        if (binNum == binRC){
            continue;
        }
        if (provisionalRepList.at(binRC) != RC){
            valid = false;
            logFile << "the reps of bins " << binNum << " and " << binRC << " are not mirrored, this points to an error" << endl;
            logFile << "rep : " << repK << " RC : " << RC << endl;
            cerr << "the reps of bins " << binNum << " and " << binRC << " are not mirrored, this points to an error" << endl;
            cerr << "rep : " << repK << " RC : " << RC << endl;

            logFile << "the Kmers in the RC bin number: " << binRC << endl;
            cerr << "the Kmers in the RC bin number: " << binRC << endl;
            for (const auto& kmer : reverseBins.at(binRC)){
                logFile << kmer << endl;
                cerr << kmer << endl;
            }
            logFile << "the Kmers in the rep bin number: " << binNum << endl;
            cerr << "the Kmers in the rep bin number: " << binNum << endl;
            for (const auto& kmer : reverseBins.at(binNum)){
                logFile << kmer << endl;
                cerr << kmer << endl;
            }
        }
        if (reverseBins.at(binNum).size() != reverseBins.at(binRC).size()){
            valid = false;
            logFile << "the expected mirror bins " << binNum << " and " << binRC << " are not the same size, this points to an error" << endl;
            logFile << "rep : " << repK << " RC : " << RC << endl;
            cerr << "the expected mirror bins " << binNum << " and " << binRC << " are not the same size, this points to an error" << endl;
            cerr << "rep : " << repK << " RC : " << RC << endl;
        }

    }
    if (valid){
        logFile << "bins are valid" << endl;
    }
}