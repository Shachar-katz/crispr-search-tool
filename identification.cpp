//
//  Step_1.cpp
//  repeatesSearchProjectAltK
//
//  Created by Shachar Katz on 12/18/24.
//

#include "pipelineSectionsHeader.h"

int identifyingRepeatPatterns(string inputFileType,
                              string outputFile,
                              int minK,
                              int minLegitimateSpacer,
                              int maxLegitimateSpacer,
                              bool strict,
                              bool preStrict,
                              int interval,
                              int maxK,
                              int segmentSize,
                              int smoothingWindow,
                              int numRepeatUnits,
                              double seedPercentage,
                              string inputFile,
                              string inputFileR2,
                              const vector<string>& inputFileList,
                              const vector<string>& inputFileListR2)
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

     // Initialize stats
    stats["number_of_reads_in_file: "] = 0;
    stats["number_of_reads_in_file_with_repeat: "] = 0;
    stats["precent_reads_in_file_with_repeat: "] = 0;

     // calculate seedK and horizon
    int seedK = minK / 2;
    int horizon = (maxK + maxLegitimateSpacer) * numRepeatUnits + maxK;

    // create file list vector 
    vector<string> filesToProcess;
    
    if (!inputFileList.empty()) {
        if (inputFileType == "fastq_dual") {
            // For dual fastq, combine R1 and R2 lists
            for (int i = 0; i < inputFileList.size(); i++) {
                filesToProcess.push_back(inputFileList[i]);
                filesToProcess.push_back(inputFileListR2[i]);
            }
        } else {
            filesToProcess = inputFileList;
        }
        logFile << "Processing " << filesToProcess.size() << " files from file list for combined catalog" << endl;
    } else {
        // Single file mode - add files to vector
        if (inputFileType == "fastq_dual") {
            filesToProcess.push_back(inputFile);
            filesToProcess.push_back(inputFileR2);
        } else {
            filesToProcess.push_back(inputFile);
        }
        logFile << "Processing single file: " << inputFile << endl;
    }

     // find kmers in the files
     for (int fileIdx = 0; fileIdx < filesToProcess.size(); fileIdx++) {
        string currentFile = filesToProcess[fileIdx];
        
        if (filesToProcess.size() > 1) {
            logFile << "Processing file " << (fileIdx + 1) << "/" << filesToProcess.size() 
                    << ": " << currentFile << endl;
            cout << "Processing file " << (fileIdx + 1) << "/" << filesToProcess.size() 
                 << ": " << currentFile << endl;
        }

        try {
            MultiFormatFileReader fileReader(currentFile, inputFileType);
            logFile << "finding Kmers in file: " << currentFile << endl;
            findKmersInFile(fileReader, globalKmerMap, seedK, minK, minLegitimateSpacer, 
                          horizon, segmentSize, smoothingWindow, stats, strict, preStrict, 
                          logFile, interval, maxK);
            logFile << "Number of unique Kmers found so far: " << globalKmerMap.size() << endl;
        } catch (const exception& e) {
            logFile << "Error processing file " << currentFile << ": " << e.what() << endl;
            cerr << "Error processing file " << currentFile << ": " << e.what() << endl;
            continue;
        }
    }
     int numReadsWithRepeats = stats["number_of_reads_in_file_with_repeat: "];
     int numReads = stats["number_of_reads_in_file: "];
     double precentReadsWithRepeat = (static_cast<double>(numReadsWithRepeats) / static_cast<double>(numReads)) * 100;
     stats["precent_reads_in_file_with_repeat: "] = precentReadsWithRepeat;


    // writing output file

    ofstream outFS1;
    outFS1.open(outputFile);
    if (!outFS1.is_open()){
         logFile << "Error: Could not open output file." << endl;
         cerr << "Error: Could not open output file." << endl;
         return -1;
    }

    logFile << "opened output file named: " << outputFile << endl;
    outFS1 << left << "repeat" << '\t' << "abundance" << endl;
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