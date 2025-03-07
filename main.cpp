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
#include "Step_1.hpp"
#include "cleaningKmers.hpp"
#include "Step_3.hpp"
#include "Params.h"
using namespace std;

int main(int argc, const char * argv[]) {
    Parameters args;

    //global args
    args.add_parser("step", new ParserInteger("Step (1, 2, or 3)"), true);
    args.add_parser("inputFileType", new ParserString("Input file type"), true);

    args.add_parser("seedK", new ParserInteger("Seed k"), true);
    args.add_parser("outputFile", new ParserFilename("Output file name"), true);

    // args maked based on step usage:
    args.add_parser("inputFile", new ParserFilename("Input file name (for single fastq)")); // for 1 & 3
    args.add_parser("inputFileR1", new ParserFilename("Input file R1 (for fastq_dual)")); // for 1 & 3
    args.add_parser("inputFileR2", new ParserFilename("Input file R2 (for fastq_dual)")); // for 1 & 3
    args.add_parser("minK", new ParserInteger("Minimum k")); // for 1 & 3
    args.add_parser("legitimateSpacer", new ParserInteger("Legitimate spacer length")); // for 1
    args.add_parser("inputFileCatalog", new ParserFilename("Input file catalog")); // for 3 & 2
    args.add_parser("alpha", new ParserInteger("Alpha (number of mutations permitted for grouping kmers)")); // for 2
    
    // Read and parse
    args.read(argc, argv);
    args.parse(); // maybe needs a true
    args.verify_mandatory();
    
    int step = args.get_int("step");

    switch(step){
        case 1:{
            if (!args.is_defined("inputFileType") ||
                !args.is_defined("outputFile") ||
                !args.is_defined("minK") ||
                !args.is_defined("seedK") ||
                !args.is_defined("legitimateSpacer")) {
                cerr << "Missing mandatory arguments for step 1." << endl;
                return 1;
            }
            string inputFileType = args.get_string("inputFileType");
            if (inputFileType == "fastq_dual") {
                if (!args.is_defined("inputFileR1") || !args.is_defined("inputFileR2")) {
                    cerr << "Missing mandatory fastq_dual input files for step 1." << endl;
                    return 1;
                }
                string inputFileR1 = args.get_string("inputFileR1");
                string inputFileR2 = args.get_string("inputFileR2");
                string outputFile = args.get_string("outputFile");
                int minK = args.get_int("minK");
                int seedK = args.get_int("seedK");
                int legitimateSpacer = args.get_int("legitimateSpacer");
                string outputFileR1 = "/Users/sarahkatz/Documents/crispr-search-tool/" + outputFile + "_R1";
                string outputFileR2 = "/Users/sarahkatz/Documents/crispr-search-tool/" + outputFile + "_R2";
                step_1(inputFileR1, inputFileType, outputFileR1, seedK, minK, legitimateSpacer);
                step_1(inputFileR2, inputFileType, outputFileR2, seedK, minK, legitimateSpacer);
            } else {
                if (!args.is_defined("inputFile")) {
                    cerr << "Missing mandatory inputFile for step 1." << endl;
                    return 1;
                }
                string inputFile = args.get_string("inputFile");
                string outputFile = args.get_string("outputFile");
                int minK = args.get_int("minK");
                int seedK = args.get_int("seedK");
                int legitimateSpacer = args.get_int("legitimateSpacer");
                outputFile = "/Users/sarahkatz/Documents/crispr-search-tool/" + outputFile;
                cout << "repeat finder initialized" << endl;
                step_1(inputFile, inputFileType, outputFile, seedK, minK, legitimateSpacer);
            }
            break;
        }
        case 2:{
            if (!args.is_defined("inputFileCatalog") ||
                !args.is_defined("outputFile") ||
                !args.is_defined("seedK") ||
                !args.is_defined("alpha")) {
                cerr << "Missing mandatory arguments for step 2." << endl;
                return 1;
            }
            string inputFileCatalog = args.get_string("inputFileCatalog");
            string outputFile = args.get_string("outputFile");
            int seedK = args.get_int("seedK");
            int alpha = args.get_int("alpha");
            cleaningKmers(inputFileCatalog, outputFile, seedK, alpha);
            break;
        }
        case 3:{
            if (!args.is_defined("inputFileType") ||
                !args.is_defined("inputFileCatalog") ||
                !args.is_defined("outputFile") ||
                !args.is_defined("seedK")) {
                cerr << "Missing mandatory arguments for step 3." << endl;
                return 1;
            }
            string inputFileType = args.get_string("inputFileType");
            // decide on dual or single fastq
            if (args.is_defined("inputFileR1") && args.is_defined("inputFileR2")) {
                string inputFileR1 = args.get_string("inputFileR1");
                string inputFileR2 = args.get_string("inputFileR2");
                string inputFileCatalog = args.get_string("inputFileCatalog");
                string outputFile = args.get_string("outputFile");
                int seedK = args.get_int("seedK");
                string outputFileR1 = "/Users/sarahkatz/Documents/crispr-search-tool/" + outputFile + "_R1";
                string outputFileR2 = "/Users/sarahkatz/Documents/crispr-search-tool/" + outputFile + "_R2";
                step_3(inputFileR1, inputFileType, inputFileCatalog, outputFileR1, seedK);
                step_3(inputFileR2, inputFileType, inputFileCatalog, outputFileR2, seedK);
            } else if (args.is_defined("inputFile")) {
                string inputFile = args.get_string("inputFile");
                string inputFileCatalog = args.get_string("inputFileCatalog");
                string outputFile = args.get_string("outputFile");
                int seedK = args.get_int("seedK");
                outputFile = "/Users/sarahkatz/Documents/crispr-search-tool/" + outputFile;
                step_3(inputFile, inputFileType, inputFileCatalog, outputFile, seedK);
            } else {
                cerr << "Missing mandatory input file for step 3" << endl;
                return 1;
            }
            break;
        }
        default:{
            cerr << "you must specify a step in the running process and provide the appropriate parameters." << endl;
            return 1;
        }
    }

    return 0;
}

