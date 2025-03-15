
#include "functions.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"

// this function populates 2 maps: 
// 1) Smap: maps from Smer to a vector of Kmers
// 2) Kmap: maps from Kmer to abundance data that is recived in the file
void catalogToSAndKMaps(ifstream& InCatalog, unordered_map<string,vector<string>>& Smap, unordered_map<string,data_t>& Kmap, int seedK, ofstream& logFile) {
    logFile << "initializing converting catalog to an Smap and Kmap" << endl; 
    string line;
    string Kmer;
    int abundance;
     // we dump the header
    getline(InCatalog, line);
    // iterate over all lines in the catalog
    while (!InCatalog.eof() && InCatalog.good()) {
        getline(InCatalog, line);
        // extract the line in a stream buffer
        istringstream iss(line);
        // extract the data (Kmer string and its abundance)
        if (iss >> Kmer >> abundance) {
            // if (!isKmer(Kmer)){
            //     continue;
            // }
            // we generate a reverse complement for the Kmer
            string reverseComp = reverseComplement(Kmer);
            // we place the Kmer and reverse complement in a map mapping from Kmer to abundance, palindrom score, and length
            Kmap[Kmer].numLines = abundance;
            Kmap[Kmer].palindromicScore = palindromicScore(Kmer);
            Kmap[Kmer].KLen = Kmer.length();
            Kmap[reverseComp].numLines = abundance;
            Kmap[reverseComp].palindromicScore = palindromicScore(reverseComp);
            Kmap[reverseComp].KLen = reverseComp.length();
            // for the length of the Kmer - the length of an Smer we extract Smers of the Kmer and its reverse complement and place them in our Smap (maps from Smer to vector of Kmers)
            for (int j = 0; j <= (Kmer.size() - seedK) ; j++){
                string Smer = Kmer.substr(j,seedK);
                string reverseSmer = reverseComp.substr(j,seedK);
                if (Smer.size() == seedK){
                    Smap[Smer].emplace_back(Kmer);
                    Smap[reverseSmer].emplace_back(reverseComp);
                }

            }
        }
    }
    logFile << "after initilizing Kmap length: " << Kmap.size() << "and Smap length: " << Smap.size() << endl; //debugging
}

// this function recives a Kmer input and generates a vector of Smers that comprize it
void findSmersVect(string shortestKmer, vector<string>& smerVect, int seedK){
    for (int j = 0; j < (shortestKmer.size() - seedK); j++){
        string SmerForSearch = shortestKmer.substr(j,seedK);
        smerVect.emplace_back(SmerForSearch);
    }
}

// this function populates the potential relations set with pairs of Kmers that are potentially related as seedK of their nucleotides are identical (they are grouped under the same Smer)
void makePotentialRelationsSet(const unordered_map<string,vector<string>>& Smap, set<pair<string,string>>& potentialRelationSet, ofstream& logFile){
    logFile << "initializing making the potential relations set" << endl; 
    // for every Smer and vector of Kmers in the Smap, if the vector contains more then 1 Kmer:
    for (const auto& [Smer, KVect] : Smap){
        if (KVect.size() > 1){
            vector<string> copyKVect = KVect;
            // we sort the Kmers from small to large
            sort(copyKVect.begin(), copyKVect.end(), [](const string &a, const string &b) {
                return (a.length() == b.length()) ? pickKey(a) < pickKey(b) : a.length() < b.length();
            });

            // for the sorted vector of Kmers we put every Kmer in a pair with all the others to signify they may be subsets of one another
            for (int i = 0; i < copyKVect.size(); i++){
                for (int j = 0; j < copyKVect.size(); j++) {
                    if (i != j){
                        potentialRelationSet.emplace(copyKVect.at(i), copyKVect.at(j));
                    }
                }
            }
        }
    }
    logFile << "after initilizing potential relations set's length: " << potentialRelationSet.size() << endl; //debugging
}

// this is a simple function that gets a pair of strings and 2 empty variables and assigns them them to the variables by length (shortest, and the other one) or if equal chooses order based on cannonic form 
void chooseShortestK (const pair<string, string>& Kpair, string& shortestK, string& otherK){
    if (Kpair.first.length() < Kpair.second.length()) {
        shortestK = Kpair.first;
        otherK = Kpair.second;
    } else if (Kpair.first.length() > Kpair.second.length()){
        shortestK = Kpair.second;
        otherK = Kpair.first;
    } else if (pickKey(Kpair.first) < pickKey(Kpair.second)){
        shortestK = Kpair.first;
        otherK = Kpair.second;
    } else{
        shortestK = Kpair.second;
        otherK = Kpair.first;
    }
}

// this function takes in 2 Kmers that are verified to be related and a bins structure and manages the binning of the Kmers together
void binRelatives(const string& shortestK, const string& otherK, DynamicBins& bins){
    // if the neither Kmers have a bin, we bin them together
    if (bins.getBin(shortestK) == -1 && bins.getBin(otherK) == -1){
        bins.autoGroup(shortestK, otherK);
    }
    // if the shortest K has a bin already, we add the other K to that bin
    else if (bins.getBin(shortestK) != -1 && bins.getBin(otherK) == -1){
        int existingBin = bins.getBin(shortestK);
        bins.addToExistingBin(otherK, existingBin);
    }
    // if the other K has a bin already, we add the shortest K to that bin
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

// This function iterates over the pairs of potentially related Kmers and verifies using the Smap and a similarity score determined by the size of the Smers and an "alotted error range" variable to determine if they are related 
// If related they are binned together, otherwise they are not.
void verifyRelation(const unordered_map<string,vector<string>>& Smap, const set<pair<string,string>>& potentialRelationSet, int seedK, DynamicBins& bins, ofstream& logFile, int alpha){
    logFile << "verifying relations" << endl;
    // for every potentially related Kmer pair:
    for (const auto& Kpair : potentialRelationSet){
        // we assign them based on length to shorter and longer
        string shortestK;
        string otherK;
        chooseShortestK(Kpair, shortestK, otherK);
        // we calculate the required number of Smers that they must NOT appear in together in the vector of Smers of the shorter Kmer in order to determine that they are NOT related
        int reqDisimilarity = alpha * seedK; // length diff between 2 k's cant be longer then alpha + 1. ??
        // we then generate a vector of all the Smers that make up the shorter K
        vector<string> smerVect;
        findSmersVect(shortestK, smerVect, seedK);
        int disimilarityScore = 0;
        bool isRelated = true;
        // then we run over this generated Smer Vector:
        for(int i = 0; i < smerVect.size(); i++){
            string currentSmer = smerVect[i];
            try{
                // we then produce for every Smer in the Smer Vector of our shortest Kmer, the vector of Kmers that it appears in on the Smap
                const auto& currKvect = Smap.at(currentSmer);
                // if we don't find the longer K in this Smer we increment the disimilarity score for this pair
                if (find(currKvect.begin(), currKvect.end(), otherK) == currKvect.end()){
                    disimilarityScore++;
                }
                // if at any point the disimilarity score passes the required disimilarity score we determine that this pair is NOT related and exit the function
                if (disimilarityScore > reqDisimilarity){
                    isRelated = false;
                    break;
                }
            }
            catch(const out_of_range& ex){
                logFile << "failed to find " << currentSmer << "in Smap" << endl;
            }
        }
        // if after the check they are still related, we bin them together using the binning function, as well as their RC's, otherwise we bin them seperately.
        if (isRelated){
            binRelatives(shortestK, otherK, bins);
            // string RCShortestK = reverseComplement(shortestK);
            // string RCOtherK = reverseComplement(otherK);
            // binRelatives(RCShortestK, RCOtherK, bins); // viable_change
        }
        else{ // maybe need to add RC's here too?
            bins.addAutoSingle(shortestK);
            bins.addAutoSingle(otherK);
        }
    }
    logFile << "grouping completed" << endl; // debugging
}

// This function iterates over the Smap and bins any Kmers that were never put into a pair and thus are not in a bin
void binSingles(const unordered_map<string,vector<string>>& Smap, DynamicBins& bins, ofstream& logFile){
    logFile << "begining binning singles" << endl; // debugging
    for (const auto& [Smer, KVect] : Smap){
        if (KVect.size() == 1){
            string Kmer = KVect[0];
            bins.addAutoSingle(Kmer);
        }
    }
    logFile << "binning singles" << endl; // debugging
}

// this function takes in 2 Kmers and the Kmap and returns the one that is of higher abundance in the Kmap
string kmerCompetition(const unordered_map<string,data_t>& Kmap, string currentRep, string auditioningKmer, ofstream& logFile){
    try {
        // we search for the abundance of both the current leader Kmer and the auditioning leader Kmer and return the one of higher abundance
        int abundanceCurrent = Kmap.at(currentRep).numLines;
        int abundanceAudition = Kmap.at(auditioningKmer).numLines;

        if (abundanceCurrent < abundanceAudition) {
            return auditioningKmer;
        }
        else if(abundanceCurrent == abundanceAudition){
            return tieBreaker(currentRep,auditioningKmer);
        }
        return currentRep;
    }
    catch (const out_of_range& ex) {
        // Check which key is missing
        if (Kmap.find(currentRep) == Kmap.end()) {
            logFile << "Key not found in Kmap: " << currentRep << endl;
            cerr << "Key not found in Kmap: " << currentRep << endl;
        }
        else if (Kmap.find(auditioningKmer) == Kmap.end()) {
            logFile << "Key not found in Kmap: " << auditioningKmer << endl;
            cerr << "Key not found in Kmap: " << auditioningKmer << endl;
        }
        else {
            // (Extremely unlikely, but handle the generic case if needed)
            logFile << "A key was not found in the map, but neither currentRep nor auditioningKmer is missing??\n";
            cerr << "A key was not found in the map, but neither currentRep nor auditioningKmer is missing??\n";
        }
        return "";
    }
}

// this function cannonicly picks a Kmer rep if their abundance is equal
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
void selectReps(unordered_map<int, string>& provisionalRepList, const unordered_map<int, vector<string>>& reverseBins, const unordered_map<string,data_t>& Kmap, ofstream& logFile){
    logFile << "begining selecting reps" << endl;
    // we iterate over the reverse bins structure that points from bin number to a vector of Kmers in that bin.
    for (const auto& [binNum, KVect] : reverseBins){
        vector <string> sortedKVect = KVect;
        // we sort the Kmers in the vector from small to large
        sort(sortedKVect.begin(), sortedKVect.end(),[](const std::string &a, const std::string &b) {
                return (a.length() < b.length()) || (a.length() == b.length() && pickKey(a) < pickKey(b));
            }
        );
        // for every Kmer:
        for (const auto& auditioningKmer : sortedKVect){
            // if that bin doesn't have a representative yet (aka index 0 of our current vector) we assign the Kmer as the representative of this bin and continue
            if (provisionalRepList.count(binNum) == 0){
                provisionalRepList[binNum] = auditioningKmer;
                continue;
            }
            // if that bin DOES have a representative and weve reached this part of the code:
            // 1) we assign the representative to a variable and check which is more abundent using Kmer competition
            // 2) if the auditioning Kmer is more abundent we assign it as the new representative, else we leave the bin as is (with the existing rep)
            string currentRep = provisionalRepList.at(binNum);
            if (kmerCompetition(Kmap, currentRep, auditioningKmer,logFile) == auditioningKmer){
                provisionalRepList[binNum] = auditioningKmer;
            }
        }
    }
    logFile << "selected final reps" << endl; // debugging
}

// a function to cannonize reps
unordered_map<int, string> reCannonization(const unordered_map<int, string>& provisionalRepList, ofstream& logFile) {
    
    unordered_map<string, int> reverseRepList;
    unordered_map<int, string> finalRepList;
    // 1) flip reps list ( Kmer -> binNum )
    for(auto& [binNum, repKmer] : provisionalRepList){
        reverseRepList[repKmer] = binNum;
    }
    // 2) iterate over the reverse list, for every Kmer create a temp var cannonized = pickkey(Kmer)
    for(auto& [repKmer, binNum] : reverseRepList){
        string cannonized = pickKey(repKmer);
        // 3) check if Kmer != cannonized && cannonized is another key in reverse rep list:
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
        else{
            logFile << "failed to cannonize the Kmer: " << repKmer << endl;
            cerr << "failed to cannonize the Kmer: " << repKmer << endl;
        }
    }
    return finalRepList;
}

// this function iterates over all of the structures needed to collect all the data and format it in an output map
void creatingOutputMap(unordered_map<string,data_t>& outputMap, unordered_map<string,data_t>& binsOutputMap, const unordered_map<int,string>& choosenReps, const unordered_map<int, vector<string>>& reverseBins, const unordered_map<string,data_t>& Kmap){
    int newBinNum = 1;
    for(const auto& [binNum, repKmer] : choosenReps){
        // filling up normal output
        auto& data = outputMap[repKmer];
        data.binNum = "K_" + to_string(newBinNum);
        data.numLines = Kmap.at(repKmer).numLines; 
        data.palindromicScore = Kmap.at(repKmer).palindromicScore;
        data.KLen = Kmap.at(repKmer).KLen;
        // filling up subKmers output
        auto& subKmersVect = reverseBins.at(binNum); 
        for (string subKmer : subKmersVect){
            auto& dataBins = binsOutputMap[subKmer];
            dataBins.binNum = "K_" + to_string(newBinNum);
            dataBins.palindromicScore = Kmap.at(subKmer).palindromicScore;
            dataBins.KLen = Kmap.at(subKmer).KLen;
        }
        // updating bin nums
        newBinNum++;
    }
}

// function that calculates palindromic score for repeat (num nucleotides thar are palindromic)
int palindromicScore(string Kmer){
    int palindromScore = 0;
    int iForward = 0;
    int iReverse = Kmer.length() - 1;
    while(iForward < iReverse){
        char nucleotide = Kmer[iForward];
        char nucleotideRe = Kmer[iReverse];
        switch(nucleotide)
        {
            case 'A':
                if (nucleotideRe == 'T'){
                    palindromScore++;
                }
                break;
            case 'G': 
                if (nucleotideRe == 'C'){
                    palindromScore++;
                }
                break;
            case 'C': 
                if (nucleotideRe == 'G'){
                    palindromScore++;
                }
                break;
            case 'T': 
                if (nucleotideRe == 'A'){
                    palindromScore++;
                }
                break;
            default: break;
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
        if (provisionalRepList.at(binRC) != RC){
            valid = false;
            logFile << "the reps of bins " << binNum << " and " << binRC << " are not mirrored, this points to an error" << endl;
            logFile << "rep : " << repK << " RC : " << RC << endl;
            cerr << "the reps of bins " << binNum << " and " << binRC << " are not mirrored, this points to an error" << endl;
            cerr << "rep : " << repK << " RC : " << RC << endl;

            logFile << "the Kmers in the RC bin number: " << binRC << endl;
            cerr << "the Kmers in the RC bin number: " << binRC << endl;
            for (const auto& Kmer : reverseBins.at(binRC)){
                logFile << Kmer << endl;
                cerr << Kmer << endl;
            }
            logFile << "the Kmers in the rep bin number: " << binNum << endl;
            cerr << "the Kmers in the rep bin number: " << binNum << endl;
            for (const auto& Kmer : reverseBins.at(binNum)){
                logFile << Kmer << endl;
                cerr << Kmer << endl;
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

// bool isKmer(string maybeKmer){
//     for (int i = 0; i < maybeKmer.size(); i++) {
//         char nucleaotid = maybeKmer[i];
//         if (nucleaotid != 'A' || nucleaotid != 'T' || nucleaotid != 'C' || nucleaotid != 'G'){
//             return false;
//         }
//     }
//     return true;
// }