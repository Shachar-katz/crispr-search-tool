//
//  Step_1.cpp
//  repeatesSearchProjectAltK
//
//  Created by Sarah Katz on 12/18/24.
//

#include "Step_1.hpp"
void step_1(string inputFile, string inputFileType, string outputFile, int seedK, int minK){
    
    unordered_map<string,int> globalKmerMap;
    
    MultiFormatFileReader fileReader(inputFile, inputFileType);
    
    cout << "finding Kmers in file" << endl;
    findKmersInFile(fileReader, globalKmerMap, seedK, minK);
    
    cout << "Number of Kmers found: " << globalKmerMap.size() << endl;
    
    ofstream outFS;
    outFS.open(outputFile);
    if (!outFS.is_open()){
         cerr << "Error: Could not open output file." << endl;
         return;
    }
    
    cout << "opened output file named: " << outputFile << endl;
    outFS << setw(30)<< left << "repeat" << '\t' << "number_of_lines" << '\n';
    writeUnorderedMapToFile(globalKmerMap, outFS);
    cout << "written" << endl;
    outFS.close();
}
