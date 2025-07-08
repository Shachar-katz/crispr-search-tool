//
//  main.cpp
//  repeatesSearchProject
//
//  Created by Sarah Katz on 10/23/24.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include "functions.hpp"
#include "classData.hpp"
#include "pipelineSectionsHeader.h"
#include "Params.h"
using namespace std;

void init_params(const char* name, int argc, const char **argv, Parameters& args)
{
      //global args
    args.add_parser("step", new ParserInteger("Step (1, 2, 3, 4 or 5(combo))"), true);
    args.add_parser("outputFile", new ParserFilename("Output file name"), true);

    // args maked based on step usage:
    args.add_parser("seedPercentage", new ParserDouble("precentage of min k that makes up the length of the seed", 0.5)); // for 1 & 2
    args.add_parser("inputFileType", new ParserString("Input file type")); // for 1 & 3
    args.add_parser("inputFile", new ParserFilename("Input file name (for single fastq)")); // for 1 & 3
    args.add_parser("inputFileList", new ParserFilename("Input file containing list of files (for combo step)")); // for 5
    args.add_parser("inputFileListR1", new ParserFilename("Input file containing list of R1 files (for fastq_dual combo)")); // for 5
    args.add_parser("inputFileListR2", new ParserFilename("Input file containing list of R2 files (for fastq_dual combo)")); // for 5    
    args.add_parser("inputFileR1", new ParserFilename("Input file R1 (for fastq_dual)")); // for 1 & 3
    args.add_parser("inputFileR2", new ParserFilename("Input file R2 (for fastq_dual)")); // for 1 & 3
    args.add_parser("minK", new ParserInteger("Minimum k", 20)); // for 1 & 3
    args.add_parser("minLegitimateSpacer", new ParserInteger("min Legitimate spacer length", 10)); // for 1 & 3
    args.add_parser("maxMismatchesForKmers", new ParserInteger("max mismatches when identifying a known kmer", 10)); // for array dump
    args.add_parser("maxLegitimateSpacer", new ParserInteger("max Legitimate spacer length", 200)); // for 1
    args.add_parser("inputFileCatalog", new ParserFilename("Input file catalog")); // for 3 and 2
    args.add_parser("secondInputFileCatalog", new ParserFilename("Additional input file catalog for dual fastq")); // for 2
    args.add_parser("alpha", new ParserInteger("Alpha (number of mutations permitted for grouping kmers)")); // for 2
    args.add_parser("interval", new ParserInteger("interval for reporting progress through the file", 0)); // for 1 & 3
    args.add_parser("maxK", new ParserInteger("maximum allowed length for a repeat", 70)); // for 1
    args.add_parser("horizonCoefficient", new ParserInteger("the coefficient to generate a search window w/i line (used to multiply the num repetative units AKA spacer-repeat sets)", 4)); // for 1
    args.add_parser("segmentSize", new ParserInteger("the size of the segments for calculating repeatition score.)", 100)); // for 1
    args.add_parser("smoothingWindow", new ParserInteger("the number of segemnts used for smoothing the repeatition score.)", 2)); // for 1
    args.add_parser("preStrict", new ParserBoolean("Should we throw suspected tandam lines before processing", false)); //for 1
    args.add_parser("strictDuring", new ParserBoolean("Should we throw suspected tandam lines during processing", false)); //for 1
    args.add_parser("repSelectWithWeight", new ParserBoolean("Should we select reps with weight (true) or with abundance", true)); //for 2
    
     if (argc == 1) {
        args.usage(name);
        return;
    }
    // Read and parse
    args.read(argc, argv);
    args.parse(); // maybe needs a true
    args.verify_mandatory();
}

bool step_1_executor(Parameters& args){
    // setting all the variables from args
    string inputFileType = args.get_string("inputFileType");
    string outputFile = args.get_string("outputFile");
    double seedPercentage = args.get_double("seedPercentage");
    int minK = args.get_int("minK");
    int minLegitimateSpacer = args.get_int("minLegitimateSpacer");
    int maxLegitimateSpacer = args.get_int("maxLegitimateSpacer");
    int interval = args.get_int("interval");
    int maxK = args.get_int("maxK");
    int numRepetativeUnits = args.get_int("horizonCoefficient");
    int segmentSize = args.get_int("segmentSize");
    int smoothingWindow = args.get_int("smoothingWindow");
    bool strict = args.get_bool("strictDuring");
    bool preStrict = args.get_bool("preStrict");
    string inputFileR1 = "";
    string inputFileR2 = "";
    // a case for dual fastq
    if (args.is_defined("inputFileR1") && args.is_defined("inputFileR2")) {
        inputFileR1 = args.get_string("inputFileR1");
        inputFileR2 = args.get_string("inputFileR2");
        if (interval == 0){ interval = 100000; }
        cout << "repeat finder initialized for R1" << endl;
    } 
    else if (args.is_defined("inputFile")) {
        inputFileR1 = args.get_string("inputFile");
        inputFileR2 = "";
        if (interval == 0){ interval = 1000; }
        cout << "repeat finder initialized" << endl;
    }
    else {
        cerr << "Missing mandatory input file/s for step 1" << endl;
        return false;
    }
    int run = identifyingRepeatPatterns(inputFileType, outputFile, minK, minLegitimateSpacer, maxLegitimateSpacer, strict, preStrict, interval, maxK, segmentSize, smoothingWindow, numRepetativeUnits, seedPercentage, inputFileR1, inputFileR2);
    if (run != 0){
        cerr << "ERROR: could not complete identifying repeat pattern run, please refer to previous error messages for more information." << endl;
        return false;
    }
    return true;
}

bool step_1_executor_combo(Parameters& args, const vector<string>& inputFiles, const string& outputFile, const vector<string>& inputFilesR2 = vector<string>()) {
    string inputFileType = args.get_string("inputFileType");
    double seedPercentage = args.get_double("seedPercentage");
    int minK = args.get_int("minK");
    int minLegitimateSpacer = args.get_int("minLegitimateSpacer");
    int maxLegitimateSpacer = args.get_int("maxLegitimateSpacer");
    int interval = args.get_int("interval");
    int maxK = args.get_int("maxK");
    int numRepetativeUnits = args.get_int("horizonCoefficient");
    int segmentSize = args.get_int("segmentSize");
    int smoothingWindow = args.get_int("smoothingWindow");
    bool strict = args.get_bool("strictDuring");
    bool preStrict = args.get_bool("preStrict");
    
    if (interval == 0) {
        interval = (inputFileType == "fastq_dual") ? 100000 : 1000;
    }
    
    // Validation logic
    if (inputFiles.empty()) {
        cerr << "Error: No input files provided" << endl;
        return false;
    }
    
    if (inputFileType == "fastq_dual") {
        if (inputFilesR2.empty() || inputFiles.size() != inputFilesR2.size()) {
            cerr << "Error: For fastq_dual, both R1 and R2 file lists must be provided and have same size" << endl;
            return false;
        }
    }
    
    cout << "Running step 1 identification on " << inputFiles.size() << " files" << endl;
    
    // Call the modified function with the file lists
    int run = identifyingRepeatPatterns(inputFileType, outputFile, minK, 
                                       minLegitimateSpacer, maxLegitimateSpacer, strict, 
                                       preStrict, interval, maxK, segmentSize, 
                                       smoothingWindow, numRepetativeUnits, seedPercentage, 
                                       "", "", inputFiles, inputFilesR2);
    
    if (run != 0) {
        cerr << "ERROR: could not complete identifying repeat pattern run for combo step" << endl;
        return false;
    }
    return true;
}

bool step_2_executor(Parameters& args){
    string inputFileCatalog;
    string outputFile;
    
    // Check if this is being called from combo step (step 5)
    if (args.get_int("step") == 5) {
        // For combo step, use the generated catalog files
        string baseOutputFile = args.get_string("outputFile");
        inputFileCatalog = baseOutputFile + "_combined_catalog";
        outputFile = baseOutputFile + "_cleaned_catalog";
    } else {
        // For regular step 2, use the provided arguments
        inputFileCatalog = args.get_string("inputFileCatalog");
        outputFile = args.get_string("outputFile");
    }
    
    string inputFileCatalog2 = "";
    if (args.is_defined("secondInputFileCatalog")) {
        inputFileCatalog2 = args.get_string("secondInputFileCatalog");
    }
    
    double seedPercentage = args.get_double("seedPercentage");
    int minK = args.get_int("minK");
    int alpha = args.get_int("alpha");
    bool weightSelector = args.get_bool("repSelectWithWeight");
    
    int run = cleaningKmers(inputFileCatalog, outputFile, minK, alpha, weightSelector, seedPercentage, inputFileCatalog2);
    if (run != 0){
        cerr << "ERROR: could not complete the catalog cleaning run, please refer to previous error messages for more information." << endl;
        return false;
    }
    return true;
}

bool step_3_executor(Parameters& args, string fileR1 = "", string fileR2 = "", string catalogFile = "", string outputFileName = ""){
    string inputFileType = args.get_string("inputFileType");
    string inputFileCatalog;
    string outputFile;
    int interval = args.get_int("interval");
    int minLegitimateSpacer = args.get_int("minLegitimateSpacer");
    int minK = args.get_int("minK");
    double seedPercentage = args.get_double("seedPercentage");
    string inputFileR1 = "";
    string inputFileR2 = "";
    
    // Check if this is being called from combo step (step 5)
    if (args.get_int("step") == 5) {
        // For combo step, use the passed parameters
        inputFileR1 = fileR1;
        inputFileR2 = fileR2;
        inputFileCatalog = catalogFile;
        outputFile = outputFileName;
        if (interval == 0){ 
            interval = (inputFileType == "fastq_dual") ? 100000 : 1000; 
        }
    } else {
        // For regular step 3, use the provided arguments
        inputFileCatalog = args.get_string("inputFileCatalog");
        outputFile = args.get_string("outputFile");
        
        // decide on dual or single fastq
        if (args.is_defined("inputFileR1") && args.is_defined("inputFileR2")) {
            inputFileR1 = args.get_string("inputFileR1");
            inputFileR2 = args.get_string("inputFileR2");
            if (interval == 0){ interval = 100000; }
            cout << "repeat finder initialized for R1" << endl;
        } 
        else if (args.is_defined("inputFile")) {
            inputFileR1 = args.get_string("inputFile");
            inputFileR2 = "";
            if (interval == 0){ interval = 1000; }
            cout << "repeat finder initialized" << endl;
        } 
        else {
            cerr << "Missing mandatory input file/s for step 3" << endl;
            return false;
        }
    }

    int run = findingKnownRepeats(inputFileR1, inputFileType, inputFileCatalog, outputFile, minLegitimateSpacer, minK, interval, seedPercentage, inputFileR2);
    if (run != 0){
        cerr << "ERROR: could not complete the finding known repeats run, please refer to previous error messages for more information." << endl;
        return false;
    }
    return true;
}

// Modified array_dump_executor
bool array_dump_executor(Parameters& args, string fileR1 = "", string fileR2 = "", string catalogFile = "", string outputFileName = ""){
    string inputFileType = args.get_string("inputFileType");
    string inputFileCatalog;
    string outputFile;
    int interval = args.get_int("interval");
    int minLegitimateSpacer = args.get_int("minLegitimateSpacer");
    int maxLegitimateSpacer = args.get_int("maxLegitimateSpacer");
    int minK = args.get_int("minK");
    int maxMismatches = args.get_int("maxMismatchesForKmers");
    double seedPercentage = args.get_double("seedPercentage");
    string inputFileR1 = "";
    string inputFileR2 = "";
    
    // Check if this is being called from combo step (step 5)
    if (args.get_int("step") == 5) {
        // For combo step, use the passed parameters
        inputFileR1 = fileR1;
        inputFileR2 = fileR2;
        inputFileCatalog = catalogFile;
        outputFile = outputFileName;
        if (interval == 0){ 
            interval = (inputFileType == "fastq_dual") ? 100000 : 1000; 
        }
    } else {
        // For regular step 4, use the provided arguments
        inputFileCatalog = args.get_string("inputFileCatalog");
        outputFile = args.get_string("outputFile");
        
        // decide on dual or single fastq
        if (args.is_defined("inputFileR1") && args.is_defined("inputFileR2")) {
            inputFileR1 = args.get_string("inputFileR1");
            inputFileR2 = args.get_string("inputFileR2");
            if (interval == 0){ interval = 100000; }
            cout << "repeat finder initialized for R1" << endl;
        } 
        else if (args.is_defined("inputFile")) {
            inputFileR1 = args.get_string("inputFile");
            inputFileR2 = "";
            if (interval == 0){ interval = 1000; }
            cout << "repeat finder initialized" << endl;
        } 
        else {
            cerr << "Missing mandatory input file for array dump" << endl;
            return false;
        }
    }
    
    int run = arrayDump(inputFileR1, inputFileType, inputFileCatalog, outputFile, minLegitimateSpacer, maxLegitimateSpacer, minK, interval, maxMismatches, seedPercentage, inputFileR2);
    if (run != 0){
        cerr << "ERROR: could not complete the array dump run, please refer to previous error messages for more information." << endl;
        return false;
    }
    return true;
}

bool combo_executor(Parameters& args){
    // Get file lists
    vector<string> inputFilesR1;
    vector<string> inputFilesR2;
    string inputFileType = args.get_string("inputFileType");
    
    if (inputFileType == "fastq_dual") {
        if (!args.is_defined("inputFileListR1") || !args.is_defined("inputFileListR2")) {
            cerr << "Missing mandatory input file lists for fastq_dual combo step" << endl;
            return false;
        }
        inputFilesR1 = readFileList(args.get_string("inputFileListR1"));
        inputFilesR2 = readFileList(args.get_string("inputFileListR2"));
        
        if (inputFilesR1.size() != inputFilesR2.size()) {
            cerr << "Error: R1 and R2 file lists must have the same number of files" << endl;
            return false;
        }
    } else {
        if (!args.is_defined("inputFileList")) {
            cerr << "Missing mandatory input file list for combo step" << endl;
            return false;
        }
        inputFilesR1 = readFileList(args.get_string("inputFileList"));
        // inputFilesR2 remains empty for non-dual fastq
    }
    
    if (inputFilesR1.empty()) {
        cerr << "Error: No input files found in file list" << endl;
        return false;
    }
    
    string baseOutputFile = args.get_string("outputFile");
    
    // STEP 1: Run identification on all files to create combined catalog
    cout << "Identifying repeat patterns across all files" << endl;
    string combinedCatalogFile = baseOutputFile + "_combined_catalog";
    
    if (!step_1_executor_combo(args, inputFilesR1, combinedCatalogFile, inputFilesR2)) {
        cerr << "ERROR: Step 1 (identification) failed" << endl;
        return false;
    }
    
    cout << "Step 1 completed. Combined catalog: " << combinedCatalogFile << endl;
    
    // STEP 2: Clean the combined catalog
    cout << "Cleaning combined catalog" << endl;
    string cleanedCatalogFile = baseOutputFile + "_cleaned_catalog";
    
    if (!step_2_executor(args)) {
        cerr << "ERROR: Step 2 (cleaning) failed" << endl;
        return false;
    }
    
    cout << "Step 2 completed. Cleaned catalog: " << cleanedCatalogFile << endl;
    
    // STEPS 3 & 4: Process each file individually
    cout << "Processing individual files" << endl;
    
    for (int i = 0; i < inputFilesR1.size(); i++) {
        string fileR1 = inputFilesR1[i];
        string fileR2 = (inputFileType == "fastq_dual") ? inputFilesR2[i] : "";
        
        // Extract base filename for output naming
        string baseFileName = fileR1.substr(fileR1.find_last_of("/\\") + 1);
        if (baseFileName.find('.') != string::npos) {
            baseFileName = baseFileName.substr(0, baseFileName.find_last_of('.'));
        }
        
        cout << "Processing file " << (i + 1) << "/" << inputFilesR1.size() << ": " << baseFileName << endl;
        
        // STEP 3: Find known repeats
        string step3OutputFile = baseOutputFile + "_" + baseFileName + "_known_repeats";
        
        if (!step_3_executor(args, fileR1, fileR2, cleanedCatalogFile, step3OutputFile)) {
            cerr << "ERROR: Step 3 failed for file: " << baseFileName << endl;
            continue; // Continue with other files
        }
        
        // STEP 4: Array dump
        string step4OutputFile = baseOutputFile + "_" + baseFileName + "_arrays";
        
        if (!array_dump_executor(args, fileR1, fileR2, cleanedCatalogFile, step4OutputFile)) {
            cerr << "ERROR: Step 4 failed for file: " << baseFileName << endl;
            continue; // Continue with other files
        }
        
        cout << "Completed processing: " << baseFileName << endl;
    }
    
    cout << "=== COMBO STEP COMPLETED ===" << endl;
    cout << "Combined catalog: " << combinedCatalogFile << endl;
    cout << "Cleaned catalog: " << cleanedCatalogFile << endl;
    cout << "Individual file outputs: " << baseOutputFile << "_[filename]_*" << endl;
    
    return true;
}

int main(int argc, const char * argv[]) {
    Parameters args;
    
    init_params(argv[0], argc, argv, args);
    
    int step = args.get_int("step");

    switch(step){
        case 1:{
            if (!args.is_defined("inputFileType") ||
                !args.is_defined("outputFile") ||
                !args.is_defined("minK") ||
                !args.is_defined("minLegitimateSpacer")) {
                cerr << "Missing mandatory arguments for step 1." << endl;
                return 1;
            }
            if(!step_1_executor(args)){ return 1; }
            break;
        }
        case 2:{
            if (!args.is_defined("inputFileCatalog") ||
                !args.is_defined("outputFile") ||
                !args.is_defined("alpha")) {
                cerr << "Missing mandatory arguments for step 2." << endl;
                return 1;
            }
            if (!step_2_executor(args)) { return 1; }
            break;
        }
        case 3:{
            if (!args.is_defined("inputFileType") ||
                !args.is_defined("inputFileCatalog") ||
                !args.is_defined("outputFile")) {
                cerr << "Missing mandatory arguments for step 3." << endl;
                return 1;
            }
            if(!step_3_executor(args)){ return 1; }
            break;
        }
        case 4:{
            if (!args.is_defined("inputFileType") ||
                !args.is_defined("inputFileCatalog") ||
                !args.is_defined("outputFile")) {
                cerr << "Missing mandatory arguments for array dump." << endl;
                return 1;
            }
            if(!array_dump_executor(args)){ return 1; }
            break;
        }
        case 5:{
            if (!args.is_defined("inputFileType") ||
                !args.is_defined("outputFile") ||
                !args.is_defined("minK") ||
                !args.is_defined("minLegitimateSpacer") ||
                !args.is_defined("alpha")) 
            {
                cerr << "Missing mandatory arguments for combo step." << endl;
                return 1;
            }
            
            string inputFileType = args.get_string("inputFileType");
            if (inputFileType == "fastq_dual") {
                if (!args.is_defined("inputFileListR1") || !args.is_defined("inputFileListR2")) {
                    cerr << "Missing inputFileListR1 and/or inputFileListR2 for fastq_dual combo step." << endl;
                    return 1;
                }
            } else {
                if (!args.is_defined("inputFileList")) {
                    cerr << "Missing inputFileList for combo step." << endl;
                    return 1;
                }
            }
            
            if(!combo_executor(args)){ return 1; }
            break;
        }
        default:{
            cerr << "you must specify a step in the running process and provide the appropriate parameters." << endl;
            return 1;
        }
    }

    return 0;
}

