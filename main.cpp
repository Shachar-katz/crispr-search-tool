//
//  main.cpp
//  repeatesSearchProject
//
//  Created by Sarah Katz on 10/23/24.
//

#include <iostream>
#include <fstream>
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
    args.add_parser("step", new ParserInteger("Step (1, 2, or 3)"), true);
    args.add_parser("outputFile", new ParserFilename("Output file name"), true);

    // args maked based on step usage:
    args.add_parser("seedPercentage", new ParserDouble("precentage of min k that makes up the length of the seed", 0.5)); // for 1 & 2
    args.add_parser("inputFileType", new ParserString("Input file type")); // for 1 & 3
    args.add_parser("inputFile", new ParserFilename("Input file name (for single fastq)")); // for 1 & 3
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
    args.add_parser("maxK", new ParserInteger("maximum allowed length for a repeat", 75)); // for 1
    args.add_parser("horizonCoefficient", new ParserInteger("the coefficient to generate a search window w/i line (used to multiply the num repetative units AKA spacer-repeat sets)", 4)); // for 1
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
    bool strict = args.get_bool("strictDuring");
    bool preStrict = args.get_bool("preStrict");
    // a case for either dual fastq
    if (inputFileType == "fastq_dual") {
        if (!args.is_defined("inputFileR1") || !args.is_defined("inputFileR2")) {
            cerr << "Missing mandatory fastq_dual input files for step 1." << endl;
            return false;
        }
        // run for strand 1
        string inputFileR1 = args.get_string("inputFileR1");
        string inputFileR2 = args.get_string("inputFileR2");
        if (interval == 0){ interval = 100000; }
        cout << "repeat finder initialized for R1" << endl;
        int run = identifyingRepeatPatterns(inputFileR1, inputFileType, outputFile, minK, minLegitimateSpacer, maxLegitimateSpacer, strict, preStrict, interval, maxK, numRepetativeUnits, seedPercentage, inputFileR2);
        if (run != 0){
            cerr << "ERROR: could not complete identifying repeat pattern run, please refer to previous error messages for more information." << endl;
            return false;
        }
    } 
    // a case for all other file types
    else {
        if (!args.is_defined("inputFile")) {
            cerr << "Missing mandatory inputFile for step 1." << endl;
            return false;
        }
        string inputFile = args.get_string("inputFile");
        if (interval == 0){ interval = 1000; }
        cout << "repeat finder initialized" << endl;
        int run = identifyingRepeatPatterns(inputFile, inputFileType, outputFile, minK, minLegitimateSpacer, maxLegitimateSpacer, strict, preStrict, interval, maxK, numRepetativeUnits, seedPercentage);
        if (run != 0){
            cerr << "ERROR: could not complete identifying repeat pattern run, please refer to previous error messages for more information." << endl;
            return false;
        }
    }
    return true;
}

bool step_2_executor(Parameters& args){
    string inputFileCatalog = args.get_string("inputFileCatalog");
    string inputFileCatalog2 = args.get_string("secondInputFileCatalog");
    string outputFile = args.get_string("outputFile");
    double seedPercentage = args.get_double("seedPercentage");
    int minK =  args.get_int("minK");
    int alpha = args.get_int("alpha");
    bool weightSelector = args.get_bool("repSelectWithWeight");
    int run = cleaningKmers(inputFileCatalog, outputFile, minK, alpha, weightSelector, seedPercentage, inputFileCatalog2);
    if (run != 0){
        cerr << "ERROR: could not complete the catalog cleaning run, please refer to previous error messages for more information." << endl;
        return false;
    }
    return true;
}

bool step_3_executor(Parameters& args){
    string inputFileType = args.get_string("inputFileType");
    string inputFileCatalog = args.get_string("inputFileCatalog");
    string outputFile = args.get_string("outputFile");
    int interval = args.get_int("interval");
    int minLegitimateSpacer = args.get_int("minLegitimateSpacer");
    int minK = args.get_int("minK");
    double seedPercentage = args.get_double("seedPercentage");
    // decide on dual or single fastq
    if (args.is_defined("inputFileR1") && args.is_defined("inputFileR2")) {
        // run for strand 1
        string inputFileR1 = args.get_string("inputFileR1");
        string inputFileR2 = args.get_string("inputFileR2");
        int run = findingKnownRepeats(inputFileR1, inputFileType, inputFileCatalog, outputFile, minLegitimateSpacer, minK, interval, seedPercentage, inputFileR2);
        if (run != 0){
            cerr << "ERROR: could not complete the finding known repeats run, please refer to previous error messages for more information." << endl;
            return false;
        }
    } 
    else if (args.is_defined("inputFile")) {
        string inputFile = args.get_string("inputFile");
        int run = findingKnownRepeats(inputFile, inputFileType, inputFileCatalog, outputFile, minLegitimateSpacer, minK, interval, seedPercentage);
        if (run != 0){
            cerr << "ERROR: could not complete the finding known repeats run, please refer to previous error messages for more information." << endl;
            return false;
        }
    } 
    else {
        cerr << "Missing mandatory input file for step 3" << endl;
        return false;
    }
    return true;
}

bool array_dump_executor(Parameters& args){
    string inputFileType = args.get_string("inputFileType");
    string inputFileCatalog = args.get_string("inputFileCatalog");
    string outputFile = args.get_string("outputFile");
    int interval = args.get_int("interval");
    int minLegitimateSpacer = args.get_int("minLegitimateSpacer");
    int maxLegitimateSpacer = args.get_int("maxLegitimateSpacer");
    int minK = args.get_int("minK");
    int maxMismatches = args.get_int("maxMismatchesForKmers");
    double seedPercentage = args.get_double("seedPercentage");
    // decide on dual or single fastq
    if (args.is_defined("inputFileR1") && args.is_defined("inputFileR2")) {
        string inputFileR1 = args.get_string("inputFileR1");
        string inputFileR2 = args.get_string("inputFileR2");
        int run = arrayDump(inputFileR1, inputFileType, inputFileCatalog, outputFile, minLegitimateSpacer, maxLegitimateSpacer, minK, interval, maxMismatches, seedPercentage, inputFileR2);
        if (run != 0){
            cerr << "ERROR: could not complete the array dump run, please refer to previous error messages for more information." << endl;
            return false;
        }
    } 
    else if (args.is_defined("inputFile")) {
        string inputFile = args.get_string("inputFile");
        int run = arrayDump(inputFile, inputFileType, inputFileCatalog, outputFile, minLegitimateSpacer, maxLegitimateSpacer, minK, interval, maxMismatches, seedPercentage);
        if (run != 0){
            cerr << "ERROR: could not complete the array dump run, please refer to previous error messages for more information." << endl;
            return false;
        }
    } 
    else {
        cerr << "Missing mandatory input file for array dump" << endl;
        return false;
    }
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
        default:{
            cerr << "you must specify a step in the running process and provide the appropriate parameters." << endl;
            return 1;
        }
    }

    return 0;
}

