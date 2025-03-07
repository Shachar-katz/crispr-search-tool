//
//  Step_3.cpp
//  repeatesSearchProjectAltK
//
//  Created by Sarah Katz on 12/29/24.
//

#include "Step_3.hpp"
#include "fileReadClass.hpp"

void step_3(string inputRead, string inputReadFileType, string inputCatalog, string outputFile, int seedK){
    unordered_map<string,Kmap_t> Smap;
    unordered_map<string,data_t> globalKmerMap;
    unordered_map<string,double> stats;
    
    ifstream catalogFile;
    catalogFile.open(inputCatalog);
    if (!isInputFileValid(catalogFile)){
        return;
    }
    cout << "catalog file opened" << endl;
    
    catalogFile.clear();
    catalogFile.seekg(0, ios::beg);
    
    buildSmap(catalogFile, Smap, seedK);
    
    cout << "smap built" << endl;
    
    catalogFile.close();
    
  //potentially add try catch for if file doesnt open
    
    MultiFormatFileReader fileReader(inputRead, inputReadFileType);
    cout << "reads file opened" << endl;
    
    findKmersInFileWithSmap(fileReader, globalKmerMap, Smap, seedK, stats);
    
    cout << "Kmers found" << endl;
    
    // writing output file
    
    ofstream outFS1;
    outFS1.open(outputFile);
    if (!outFS1.is_open()){
         cerr << "Error: Could not open output file." << endl;
         return;
    }
    
    cout << "opened output file named: " << outputFile << endl;
    outFS1 << "repeat" << '\t' << "repeats_in_file" << '\t' << "number_of_lines" << '\n';
    writeUnorderedMapToFile(globalKmerMap, outFS1);
    cout << "written" << endl;
    outFS1.close();

        // writing statistics file

    ofstream outFS2;
    string statsOutput = "stats_step_3";
    outFS2.open(statsOutput);
    if (!outFS2.is_open()){
         cerr << "Error: Could not open stats output file." << endl;
         return;
    }
    
    cout << "opened output file named: " << statsOutput << endl;
    writeUnorderedMapToFile(stats, outFS2);
    cout << "written" << endl;
    outFS2.close();
}
