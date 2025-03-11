//
//  Step_1.cpp
//  repeatesSearchProjectAltK
//
//  Created by Sarah Katz on 12/18/24.
//

#include "Step_1.hpp"
void step_1(string inputFile, string inputFileType, string outputFile, int seedK, int minK, 
            int legitimateSpacer, bool strict, bool preStrict){
    // open log file
    ofstream logFile;
    string logFileName = outputFile + "_run_log";
    logFile.open(logFileName);
    if (!logFile.is_open()){
         cerr << "Error: Could not open log output file." << endl;
         return;
    }

    unordered_map<string,int> globalKmerMap;
    unordered_map<string,double> stats;
    MultiFormatFileReader fileReader(inputFile, inputFileType);
    
    logFile << "finding Kmers in file" << endl;
    findKmersInFile(fileReader, globalKmerMap, seedK, minK, legitimateSpacer, stats, strict, preStrict, logFile);
    logFile << "Number of Kmers found: " << globalKmerMap.size() << endl;
    
    // writing output file 

    ofstream outFS1;
    outFS1.open(outputFile);
    if (!outFS1.is_open()){
         logFile << "Error: Could not open output file." << endl;
         cerr << "Error: Could not open output file." << endl;
         return;
    }
    
    logFile << "opened output file named: " << outputFile << endl;
    outFS1 << setw(30)<< left << "repeat" << '\t' << "number_of_lines" << '\n';
    writeUnorderedMapToFile(globalKmerMap, outFS1);
    logFile << "written" << endl;
    outFS1.close();

    // writing statistics file

    ofstream outFS2;
    string statsOutput = outputFile + "_stats_step_1";
    outFS2.open(statsOutput);
    if (!outFS2.is_open()){
         logFile << "Error: Could not open stats output file." << endl;
         cerr << "Error: Could not open stats output file." << endl;
         return;
    }
    
    logFile << "opened output file named: " << statsOutput << endl;
    writeUnorderedMapToFile(stats, outFS2);
    logFile << "written" << endl;
    outFS2.close();
}
