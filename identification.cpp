//
//  Step_1.cpp
//  repeatesSearchProjectAltK
//
//  Created by Shachar Katz on 12/18/24.
//

#include "pipelineSectionsHeader.h"

int identifyingRepeatPatterns(string inputFile,
                              string inputFileType,
                              string outputFile,
                              int minK,
                              int minLegitimateSpacer,
                              int maxLegitimateSpacer,
                              bool strict,
                              bool preStrict,
                              int interval,
                              int maxK,
                              int numRepeatUnits,
                              string inputFileR2)
{
    // open log file
    ofstream logFile;
    string logFileName = outputFile + "_run_log";
    logFile.open(logFileName);
    if (!logFile.is_open()){
         cerr << "Error: Could not open log output file." << endl;
         return -1;
    }

    unordered_map<string,int> globalKmerMap;
    unordered_map<string,double> stats;
    MultiFormatFileReader fileReader(inputFile, inputFileType);

    int seedK = minK / 2;
    int horizon = (maxK + maxLegitimateSpacer) * numRepeatUnits + maxK;

    logFile << "finding Kmers in file" << endl;
    findKmersInFile(fileReader, globalKmerMap, seedK, minK, minLegitimateSpacer, maxLegitimateSpacer, horizon, stats, strict, preStrict, logFile, interval, maxK);
    logFile << "Number of Kmers found: " << globalKmerMap.size() << endl;

    if (inputFileType == "fastq_dual"){
          MultiFormatFileReader fileReader2(inputFileR2, inputFileType);
          logFile << "finding Kmers in file R2" << endl;
          findKmersInFile(fileReader2, globalKmerMap, seedK, minK, minLegitimateSpacer, maxLegitimateSpacer, horizon, stats, strict, preStrict, logFile, interval, maxK);
          logFile << "Number of Kmers found after second round: " << globalKmerMap.size() << endl;
    }
    // writing output file

    ofstream outFS1;
    outFS1.open(outputFile);
    if (!outFS1.is_open()){
         logFile << "Error: Could not open output file." << endl;
         cerr << "Error: Could not open output file." << endl;
         return -1;
    }

    logFile << "opened output file named: " << outputFile << endl;
    outFS1 << left << "repeat" << '\t' << "number_of_lines" << endl;
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
         return -1;
    }

    logFile << "opened output file named: " << statsOutput << endl;
    writeUnorderedMapToFile(stats, outFS2);
    logFile << "written" << endl;
    outFS2.close();
    return 0;
}