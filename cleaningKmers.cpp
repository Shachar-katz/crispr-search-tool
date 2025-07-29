//
//  cleaningKmers.cpp
//  repeatesSearchProjectAltK
//
//  Created by Sarah Katz on 1/2/25.
//

#include "pipelineSectionsHeader.h"

int cleaningKmers(string inputCatalog, 
                  string outputFile,
                  int minK, 
                  int alpha,
                  bool weightSelector,
                  double seedPercentage,
                  string inputCatalog2)
{
    // open log file
    ofstream logFile;
    string logFileName = outputFile + "_run_log";
    logFile.open(logFileName);
    if (!logFile.is_open()){
         cerr << "Error: Could not open log output file." << endl;
         return -1;
    }

    // formatting data from catalog
    unordered_map<string,vector<string>> smap;
    unordered_map<string,data_t> kmap;
    
    ifstream catalogFile;
    catalogFile.open(inputCatalog);
    if (!isInputFileValid(catalogFile, inputCatalog)){ return -1; }

    int seedK = minK * seedPercentage;
        
    catalogToSAndKMaps(catalogFile, smap, kmap, seedK, logFile);
    
    catalogFile.close();

    if (inputCatalog2 != ""){
        ifstream catalogFile2;
        catalogFile2.open(inputCatalog2);
        if (!isInputFileValid(catalogFile2, inputCatalog2)){ return -1; }
        
        catalogToSAndKMaps(catalogFile2, smap, kmap, seedK, logFile);
        
        catalogFile2.close();
    }

    // preforming binning and selecting most common variation of repeat
    set<pair<string,string>> potentialRelationSet;
    DynamicBins bins(1000, seedK);

    makePotentialRelationsSet(smap, potentialRelationSet, logFile);
    verifyRelation(smap, potentialRelationSet, seedK, bins, logFile, alpha);
    logFile << "size of bins before binning singles: " << bins.getLen() << endl; // debugg
    binSingles(smap, bins, logFile);
    logFile << "size of bins: " << bins.getLen() << endl; // debugg
    bins.forceReclusterAll();
    bins.reverseBin();

    unordered_map<int,vector<string>> reverseBins = bins.getReBins();
    unordered_map<int,string> provisionalReps;
    
    if (weightSelector == true){
        selectRepsWeight(provisionalReps, reverseBins, kmap, seedK, logFile);
    }
    else{
        selectReps(provisionalReps, reverseBins, kmap, logFile);
    }
    
    logFile << "size of choosen reps before cannonization: " << provisionalReps.size() << endl; // debugg
    try{
        validateBins(provisionalReps,bins,reverseBins,logFile);
    }
    catch(const out_of_range& ex){
        logFile << "error occured validating bins for this data set" << endl;
        cerr << "error occured validating bins for this data set" << endl;
        return -1;
    }
    // db
    ofstream dumpBinsMissmatch;
    string dumpBinsMissmatchName = outputFile + "_error_dump_bins";
    dumpBinsMissmatch.open(dumpBinsMissmatchName);
    if (!dumpBinsMissmatch.is_open()){
         cerr << "Error: Could not open log output file." << endl;
         return -1;
    }

    unordered_map<string,int> binsMap;
    int nxtBinNum = bins.getNextBinNum();
    reSort(provisionalReps, binsMap, dumpBinsMissmatch, seedK, (alpha * seedK), nxtBinNum);
    dumpBinsMissmatch.close();
    // db
    unordered_map<int,string> finalReps = reCannonization(provisionalReps,bins,logFile);
    logFile << "size of choosen reps after cannonization: " << finalReps.size() << endl; // debugg

    // formatting output
    unordered_map<string,data_t> outputMap;
    unordered_map<string,data_t> binsOutputMap;
    string binsOutputFile = outputFile + "_bins_data";
    creatingOutputMap(outputMap, binsOutputMap, finalReps, reverseBins, kmap);
    
    ofstream outFS1;
    ofstream outFS2;
    // cout << " -" << outputFile << "- " << endl; //db
    // cout << " -" << binsOutputFile << "- " << endl; //db
    outFS1.open(outputFile);
    outFS2.open(binsOutputFile);
    if (!outFS1.is_open()) {
        logFile << "Error: Could not open output file: " << outputFile << endl;
        cerr << "Error: Could not open output file: " << outputFile << endl;
        return -1;
    }
    if (!outFS2.is_open()) {
        logFile << "Error: Could not open bins output file: " << binsOutputFile << endl;
        cerr << "Error: Could not open bins output file: " << binsOutputFile << endl;
        return -1;
    }

    logFile << "opened two output file named: " << endl << outputFile << endl << binsOutputFile << endl;
    
    outFS1 << left << "repeat" <<
     '\t' << "number_of_lines" <<
     '\t' << "bin_identifier" <<
     '\t' << "num_palindromic_nucleotides" <<
     '\t' << "kmer_Length" << endl;
    writeUnorderedMapToFile(outputMap, outFS1);
    outFS1.close();

    outFS2 << 
     left << "repeat" <<
     '\t' << "number_of_lines" <<
     '\t' << "bin_identifier" <<
     '\t' << "num_palindromic_nucleotides" <<
     '\t' << "kmer_Length" << endl;
    writeUnorderedMapToFile(binsOutputMap, outFS2);
    outFS2.close();
    return 0;
}
