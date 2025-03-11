//
//  cleaningKmers.cpp
//  repeatesSearchProjectAltK
//
//  Created by Sarah Katz on 1/2/25.
//

#include "cleaningKmers.hpp"

void cleaningKmers(string inputCatalog, string outputFile, int seedK, int alpha){
    // formatting data from catalog
    unordered_map<string,vector<string>> Smap;
    unordered_map<string,data_t> Kmap;
    
    ifstream catalogFile;
    catalogFile.open(inputCatalog);
    if (!isInputFileValid(catalogFile)){
        return;
    }
    
    catalogFile.clear();
    catalogFile.seekg(0, ios::beg);
    
    catalogToSAndKMaps(catalogFile, Smap, Kmap, seedK);
    
    catalogFile.close();

    // preforming binning and selecting most common variation of repeat
    set<pair<string,string>> potentialRelationSet;
    DynamicBins bins;
    
    makePotentialRelationsSet(Smap, potentialRelationSet);
    verifyRelation(Smap, potentialRelationSet, seedK, bins, alpha);
    cout << "size of bins before binning singles: " << bins.getLen() << endl; // debugg
    binSingles(Smap, bins);
    cout << "size of bins: " << bins.getLen() << endl; // debugg
    bins.reverseBin();

    unordered_map<int,vector<string>> reverseBins = bins.getReBins();
    unordered_map<int,string> provisionalReps;
    
    selectReps(provisionalReps, reverseBins, Kmap);
    cout << "size of choosen reps before cannonization: " << provisionalReps.size() << endl; // debugg
    try{
        validateBins(provisionalReps,bins,reverseBins);
    }
    catch(const out_of_range& ex){
        cerr << "error occured validating bins for this data set" << endl;
        return;
    }
    unordered_map<int,string> finalReps = reCannonization(provisionalReps);
    cout << "size of choosen reps after cannonization: " << finalReps.size() << endl; // debugg

    // formatting output
    unordered_map<string,data_t> outputMap;
    unordered_map<string,data_t> binsOutputMap;
    string binsOutputFile = outputFile + "_bins_data";
    creatingOutputMap(outputMap, binsOutputMap, finalReps, reverseBins, Kmap);
    
    ofstream outFS1;
    ofstream outFS2;
    // cout << " -" << outputFile << "- " << endl; //db
    // cout << " -" << binsOutputFile << "- " << endl; //db
    outFS1.open(outputFile);
    outFS2.open(binsOutputFile);
    if (!outFS1.is_open()) {
        cerr << "Error: Could not open output file: " << outputFile << endl;
    }
    if (!outFS2.is_open()) {
        cerr << "Error: Could not open bins output file: " << binsOutputFile << endl;
    }

    cout << "opened two output file named: " << endl << outputFile << endl << binsOutputFile << endl;
    
    outFS1 << 
    setw(10) << left << "repeat" <<
    setw(10) << '\t' << "number_of_lines" <<
    setw(10) << '\t' << "bin_#" <<
    setw(10) << '\t' << "num_palindromic_nucleotides" <<
    setw(10) << '\t' << "Kmer_Length" << endl;
    writeUnorderedMapToFile(outputMap, outFS1);
    outFS1.close();

    outFS2 << 
    setw(10) << left << "repeat" <<
    setw(10) << '\t' << "bin_#" <<
    setw(10) << '\t' << "num_palindromic_nucleotides" <<
    setw(10) << '\t' << "Kmer_Length" << endl;
    writeUnorderedMapToFile(binsOutputMap, outFS2);
    outFS2.close();
}
