#ifndef pipelineSectionsHeader_h
#define pipelineSectionsHeader_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include "functions.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"
#include "fileReadClass.hpp"
#include "globalFunctions.h"

using namespace std;


void step_1(string inputFile, string inputFileType, string outputFile, int seedK, int minK, int legitimateSpacer, 
            bool strict = false, bool preStrict = false, string inputFileR2 = "");

void cleaningKmers(string inputCatalog, string outputFile, int seedK, int alpha, string inputCatalog2 = "");

void step_3(string inputRead, string inputReadFileType, string inputCatalog, string outputFile, int seedK, string inputFileR2 = "");

#endif /* pipelineSectionsHeader_h */