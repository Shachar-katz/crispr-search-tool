#include "pipelineSectionsHeader.h"
#include "fileReadClass.hpp"
#include "functions.hpp"

int spacerDetection(string inputRead, 
                    string inputReadFileType, 
                    string inputCatalog, 
                    string outputFile, 
                    int minLegitimateSpacer,
                    int maxLegitimateSpacer,
                    int minK,
                    int interval,
                    int maxMismatches,
                    double seedPercentage,
                    string inputFileR2)
{
    // generating a seedK
    int seedK = minK * seedPercentage;
    int minPalindromic = 4; // turn into arg
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

    // building Smap and Kmer -> id map
    ifstream catalogFile;
    catalogFile.open(inputCatalog);
    if (!isInputFileValid(catalogFile, inputCatalog)){ return -1; }
    logFile << "catalog file opened" << endl;

    unordered_map<string,Kmap_t> smap;
    int validSmap = buildSmap(catalogFile, smap, seedK, minPalindromic);

    if (validSmap == 1){
        logFile << "header error: no header or incorrect header in the catalog file provided" << endl;
        cerr << "header error: no header or incorrect header in the catalog file provided" << endl;
        return -1;
    }
    
    logFile << "smap built" << endl;

    catalogFile.clear();
    catalogFile.seekg(0, std::ios::beg);

    unordered_map<string,string> kmerToId;
    int validKmap = buildKmap(catalogFile, kmerToId, minPalindromic);

    catalogFile.close();

    // setting up data structure for the arrays and the file reader    
    unordered_map<string,Array> globalArrayMap;
    MultiFormatFileReader fileReaderR1(inputRead, inputReadFileType);
    logFile << "reads file opened" << endl;

    // run the array identifior function
    arrayIdentifior(fileReaderR1, globalArrayMap, smap, kmerToId, seedK, stats, logFile, minLegitimateSpacer, maxLegitimateSpacer, minK, interval, maxMismatches);
    logFile << globalArrayMap.size() << "Arrays found" << endl;

    if (inputReadFileType == "fastq_dual"){
        MultiFormatFileReader fileReaderR2(inputFileR2, inputReadFileType);
        logFile << "reads file R2 opened" << endl;
        arrayIdentifior(fileReaderR2, globalArrayMap, smap, kmerToId, seedK, stats, logFile, minLegitimateSpacer, maxLegitimateSpacer, minK, interval, maxMismatches);
        logFile << globalArrayMap.size() << "Arrays found" << endl;
    }
    
    // writing output file 
    
    ofstream outFS1;
    string tableFile = outputFile + "_table";
    outFS1.open(tableFile);
    if (!outFS1.is_open()){
         logFile << "Error: Could not open table output file." << endl;
         cerr << "Error: Could not open table output file." << endl;
         return -1;
    }

    logFile << "opened output file named: " << tableFile << endl;
    outFS1 << "array_id" << '\t' << "repeat_id" << '\t' << "spacer_count" << '\t' << "array_len" << endl;;
    writeUnorderedMapToFile(globalArrayMap, outFS1);
    logFile << "written" << endl;
    outFS1.close();

    ofstream outFS2;
    string viewFile = outputFile + "_view";
    outFS2.open(viewFile);
    if (!outFS2.is_open()){
         logFile << "Error: Could not open view output file." << endl;
         cerr << "Error: Could not open view output file." << endl;
         return -1;
    }

    logFile << "opened output file named: " << viewFile << endl;
    for (auto& [id, array] : globalArrayMap){
        outFS2 << ">" << id << endl;
        outFS2 << array.getArrayStr() << endl;
    }
    logFile << "written" << endl;
    outFS2.close();

    ofstream outFS3;
    string tableFile2 = outputFile + "_table_error";
    outFS3.open(tableFile2);
    if (!outFS3.is_open()){
         logFile << "Error: Could not open table output 2 file." << endl;
         cerr << "Error: Could not open table output 2 file." << endl;
         return -1;
    }

    unordered_map<string,RepeatData> repeatMap;
    constructRepeatMap(globalArrayMap, repeatMap);

    logFile << "opened output file named: " << tableFile << endl;
    outFS3 << "array_id" << '\t' << "repeat_id" << '\t' << "repeatIdx" << '\t' << "nun_mismatches" << '\t' << "repeat_len" << endl;;
    writeUnorderedMapToFile(repeatMap, outFS3);
    logFile << "written" << endl;
    outFS3.close();

    // writing statistics file

    ofstream outFS4;
    string statsOutput = outputFile + "_stats_step_4";
    outFS4.open(statsOutput);
    if (!outFS4.is_open()){
         logFile << "Error: Could not open stats output file." << endl;
         cerr << "Error: Could not open stats output file." << endl;
         return -1;
    }
    
    logFile << "opened output file named: " << statsOutput << endl;
    writeUnorderedMapToFile(stats, outFS4);
    logFile << "written" << endl;
    outFS4.close();
    return 0;
}
