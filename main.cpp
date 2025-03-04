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
using namespace std;

int main(int argc, const char * argv[]) {
    // cout << argv[0] << endl; //(to get executable path)
    int step = atoi(argv[1]);
    switch(step){
        case 1:{
            string inputFileType = argv[2];
            if (argc < 7){
                cerr << "Error: missing args for step 1" << endl;
                cerr << "enter arguments in this order: \nexecutable, step, input_file_type, input_file_name/s, output_file_name, minimum_k, seed_k " << endl;
                return 1;
            }
            if(argc == 8 && inputFileType == "fastq_dual"){
                string inputFileR1 = argv[3];
                string inputFileR2 = argv[4];
                string outputFile = argv[5];
                int minK = atoi(argv[6]);
                int seedK = atoi(argv[7]);
                string outputFileR1 = "/Users/sarahkatz/Documents/repeat-search/repeatesSearchProjectAltK/" + outputFile + "_R1";
                string outputFileR2 = "/Users/sarahkatz/Documents/repeat-search/repeatesSearchProjectAltK/" + outputFile + "_R2";
                step_1(inputFileR1, inputFileType, outputFileR1, seedK, minK);
                step_1(inputFileR1, inputFileType, outputFileR2, seedK, minK);
            }
            else{
                string inputFile = argv[3];
                string outputFile = argv[4];
                int minK = atoi(argv[5]);
                int seedK = atoi(argv[6]);
                outputFile = "/Users/sarahkatz/Documents/repeat-search/repeatesSearchProjectAltK/" + outputFile;
                cout << "repeat finder initilized" << endl;
                step_1(inputFile, inputFileType, outputFile, seedK, minK);
            }
            break;
        }
        case 2:{
            string inputFileCatalog = argv[2];
            string outputFile2 = argv[3];
            int seedK = atoi(argv[4]);
            int alpha = atoi(argv[5]);
            if (argc < 6 || alpha < 0 || seedK < 1){
                cerr << "Error: missing or incorrect arguments for step 2" << endl;
                cerr << "enter arguments in this order: \nexecutable, step, input_file_catalog_name, output_file_name, seedK, alpha(num_mutations_premited_for_grouping_kmers)" << endl;
                return 1;
            }
            cleaningKmers(inputFileCatalog, outputFile2, seedK, alpha);
            break;
        }
        case 3:{
            string inputFileType3 = argv[2];
            if (argc < 7){
                cerr << "Error: missing args for step 3" << endl;
                cerr << "enter arguments in this order: \nexecutable, step, input_file_type, input_file_name/s, input_file_catalog_name, output_file_name, seed_k " << endl;
                return 1;
            }
            if(argc == 8 && inputFileType3 == "fastq_dual"){
                string inputFile3R1 = argv[3];
                string inputFile3R2 = argv[4];
                string inputFileCatalog3 = argv[5];
                string outputFile3 = argv[6];
                int seedK3 = atoi(argv[7]);
                string outputFileR1 = "/Users/sarahkatz/Documents/repeat-search/repeatesSearchProjectAltK/" + outputFile3 + "_R1";
                string outputFileR2 = "/Users/sarahkatz/Documents/repeat-search/repeatesSearchProjectAltK/" + outputFile3 + "_R2";
                step_3(inputFile3R1, inputFileType3, inputFileCatalog3, outputFileR1, seedK3);
                step_3(inputFile3R2, inputFileType3, inputFileCatalog3, outputFileR2, seedK3);
            }
            else{
                string inputFile3 = argv[3];
                string inputFileCatalog3 = argv[4];
                string outputFile3 = argv[5];
                int seedK3 = atoi(argv[6]);
                outputFile3 = "/Users/sarahkatz/Documents/repeat-search/repeatesSearchProjectAltK/" + outputFile3;
                step_3(inputFile3, inputFileType3, inputFileCatalog3, outputFile3, seedK3);
            }
            break;
        }
        default:
            cerr << "you must specify a step in the running process and provide the appropriate parameters." << endl;
            break;
    }

    return 0;
}

