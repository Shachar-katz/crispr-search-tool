//
//  Step_2.hpp
//  repeatesSearchProjectAltK
//
//  Created by Sarah Katz on 1/2/25.
//

#ifndef cleaningKmers_hpp
#define cleaningKmers_hpp

#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include "functions.hpp"
#include "classData.hpp"
using namespace std;

void cleaningKmers(string inputCatalog, string outputFile, int seedK, int alpha, string inputCatalog2 = "");
#endif /* cleaningKmers_hpp */
