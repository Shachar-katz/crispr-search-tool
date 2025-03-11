//
//  Step_1.cpp
//  repeatesSearchProjectAltK
//
//  Created by Sarah Katz on 12/18/24.
//

#include "Step_1.hpp"
void step_1(string inputFile, string inputFileType, string outputFile, int seedK, int minK, 
            int legitimateSpacer, bool strict, bool preStrict){
    
    unordered_map<string,int> globalKmerMap;
    unordered_map<string,double> stats;
    
    MultiFormatFileReader fileReader(inputFile, inputFileType);
    
    cout << "finding Kmers in file" << endl;
    findKmersInFile(fileReader, globalKmerMap, seedK, minK, legitimateSpacer, stats, strict, preStrict);
    
    cout << "Number of Kmers found: " << globalKmerMap.size() << endl;
    
    // writing output file 

    ofstream outFS1;
    outFS1.open(outputFile);
    if (!outFS1.is_open()){
         cerr << "Error: Could not open output file." << endl;
         return;
    }
    
    cout << "opened output file named: " << outputFile << endl;
    outFS1 << setw(30)<< left << "repeat" << '\t' << "number_of_lines" << '\n';
    writeUnorderedMapToFile(globalKmerMap, outFS1);
    cout << "written" << endl;
    outFS1.close();

    // writing statistics file

    ofstream outFS2;
    string statsOutput = "/Users/sarahkatz/Documents/data/stats_step_1";
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
