//
//  Step_3.cpp
//  repeatesSearchProjectAltK
//
//  Created by Sarah Katz on 12/29/24.
//

#include "pipelineSectionsHeader.h"
#include "fileReadClass.hpp"

int findingKnownRepeats(string inputRead, 
            string inputReadFileType, 
            string inputCatalog, 
            string outputFile, 
            int seedK, 
            int legitimateSpacer,
            int minK,
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

    // initilize stats
    unordered_map<string,double> stats;
    stats["number_of_reads_in_file: "] = 0;
    stats["number_of_reads_in_file_with_repeat: "] = 0;
    stats["precent_reads_in_file_with_repeat: "] = 0;

    
    ifstream catalogFile;
    catalogFile.open(inputCatalog);
    if (!isInputFileValid(catalogFile, inputCatalog)){ return -1; }
    logFile << "catalog file opened" << endl;

    unordered_map<string,Kmap_t> smap;
    int validSmap = buildSmap(catalogFile, smap, seedK);

    if (validSmap == 1){
        logFile << "header error: no header or incorrect header in the catalog file provided" << endl;
        cerr << "header error: no header or incorrect header in the catalog file provided" << endl;
        return -1;
    }
    
    logFile << "smap built" << endl;
    catalogFile.close();
    
  //potentially add try catch for if file doesnt open //here
    
    unordered_map<string,data_t> globalKmerMap;
    MultiFormatFileReader fileReaderR1(inputRead, inputReadFileType);
    logFile << "reads file opened" << endl;
    findKmersInFileWithSmap(fileReaderR1, globalKmerMap, smap, seedK, stats, logFile, legitimateSpacer, minK);
    logFile << globalKmerMap.size() << "Kmers found" << endl;

    if (inputReadFileType == "fastq_dual"){
        MultiFormatFileReader fileReaderR2(inputFileR2, inputReadFileType);
        logFile << "reads file R2 opened" << endl;
        findKmersInFileWithSmap(fileReaderR2, globalKmerMap, smap, seedK, stats, logFile, legitimateSpacer, minK);
        logFile << globalKmerMap.size() << "Kmers found" << endl;
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
    outFS1 << "repeat" << '\t' << "repeats_in_file" << '\t' << "number_of_lines" << endl;;
    writeUnorderedMapToFile(globalKmerMap, outFS1);
    logFile << "written" << endl;
    outFS1.close();

        // writing statistics file

    ofstream outFS2;
    string statsOutput = outputFile + "_stats_step_3";
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
