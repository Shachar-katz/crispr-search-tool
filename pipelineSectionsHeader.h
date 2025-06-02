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
#include "arrayClass.h"

using namespace std;

int identifyingRepeatPatterns(string inputFile,
                              string inputFileType,
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
                              double seedPercentage = 0.5,
                              string inputFileR2 = "");

int cleaningKmers(string inputCatalog,
                  string outputFile,
                  int minK,
                  int alpha,
                  bool weightSelector,
                  double seedPercentage = 0.5,
                  string inputCatalog2 = "");

int findingKnownRepeats(string inputRead,
                        string inputReadFileType,
                        string inputCatalog,
                        string outputFile,
                        int minLegitimateSpacer,
                        int minK,
                        int interval,
                        double seedPercentage = 0.5,
                        string inputFileR2 = "");

int arrayDump(string inputRead,
              string inputReadFileType,
              string inputCatalog,
              string outputFile,
              int minLegitimateSpacer,
              int maxLegitimateSpacer,
              int minK,
              int interval,
              int maxMismatches,
              double seedPercentage = 0.5,
              string inputFileR2 = "");

#endif /* pipelineSectionsHeader_h */