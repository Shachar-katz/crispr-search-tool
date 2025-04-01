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


int identifyingRepeatPatterns(string inputFile, string inputFileType, string outputFile, int seedK, int minK, int legitimateSpacer, 
            bool strict = false, bool preStrict = false, string inputFileR2 = "");

int cleaningKmers(string inputCatalog, string outputFile, int seedK, int alpha, string inputCatalog2 = "");

int findingKnownRepeats(string inputRead, string inputReadFileType, string inputCatalog, string outputFile, int seedK, string inputFileR2 = "", int legitimateSpacer = 10, int minK = 20);

#endif /* pipelineSectionsHeader_h */