#include "pipelineSectionsHeader.h"
#include "fileReadClass.hpp"
#include "functions.hpp"

int arrayDump(string inputRead, 
            string inputReadFileType, 
            string inputCatalog, 
            string outputFile, 
            int seedK, 
            int minLegitimateSpacer,
            int maxLegitimateSpacer,
            int minK,
            int interval,
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
    stats["number_of_reads_in_file_with_arrays: "] = 0;
    stats["precent_reads_in_file_with_arrays: "] = 0;

    
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

    catalogFile.clear();
    catalogFile.seekg(0, std::ios::beg);

    unordered_map<string,string> kmerToId;
    int validKmap = buildKmap(catalogFile, kmerToId);

    catalogFile.close();
    
  //potentially add try catch for if file doesnt open
    
    unordered_map<string,Array> globalArrayMap;
    MultiFormatFileReader fileReaderR1(inputRead, inputReadFileType);
    logFile << "reads file opened" << endl;
    arrayIdentifior(fileReaderR1, globalArrayMap, smap, kmerToId, seedK, stats, logFile, minLegitimateSpacer, maxLegitimateSpacer, minK, interval);
    logFile << globalArrayMap.size() << "Kmers found" << endl;

    if (inputReadFileType == "fastq_dual"){
        MultiFormatFileReader fileReaderR2(inputFileR2, inputReadFileType);
        logFile << "reads file R2 opened" << endl;
        arrayIdentifior(fileReaderR2, globalArrayMap, smap, kmerToId, seedK, stats, logFile, minLegitimateSpacer, maxLegitimateSpacer, minK, interval);
        logFile << globalArrayMap.size() << "Kmers found" << endl;
    }
    
    // writing output file 
    
    ofstream outFS1;
    string tableFile = outputFile + "_table";
    outFS1.open(tableFile);
    if (!outFS1.is_open()){
         logFile << "Error: Could not open output file." << endl;
         cerr << "Error: Could not open output file." << endl;
         return -1;
    }

    logFile << "opened output file named: " << tableFile << endl;
    outFS1 << "array_id" << '\t' << "repeat_id" << '\t' << "spacer_count" << endl;;
    writeUnorderedMapToFile(globalArrayMap, outFS1);
    logFile << "written" << endl;
    outFS1.close();

    ofstream outFS2;
    string viewFile = outputFile + "_view";
    outFS2.open(viewFile);
    if (!outFS2.is_open()){
         logFile << "Error: Could not open output file." << endl;
         cerr << "Error: Could not open output file." << endl;
         return -1;
    }

    logFile << "opened output file named: " << viewFile << endl;
    for (auto& [id, array] : globalArrayMap){
        outFS2 << ">" << id << endl;
        outFS2 << array.getArrayStr() << endl;
    }
    logFile << "written" << endl;
    outFS2.close();

    // writing statistics file @to here

    ofstream outFS3;
    string statsOutput = outputFile + "_stats_step_3";
    outFS3.open(statsOutput);
    if (!outFS3.is_open()){
         logFile << "Error: Could not open stats output file." << endl;
         cerr << "Error: Could not open stats output file." << endl;
         return -1;
    }
    
    logFile << "opened output file named: " << statsOutput << endl;
    writeUnorderedMapToFile(stats, outFS3);
    logFile << "written" << endl;
    outFS3.close();
    return 0;
}
