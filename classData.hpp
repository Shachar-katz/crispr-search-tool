//
//  classData.hpp
//  repeatesSearchProject
//
//  Created by Sarah Katz on 10/30/24.
//

#ifndef classData_hpp
#define classData_hpp

#include <stdio.h>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <set> 
using namespace std;

// class data to contain various data 
class data_t {
public:
    //int index = -1;
    int countInLine = 0;
    int countInFile = 0;
    int numLines = 0;
    int palindromicScore = 0;
    int kLen = 0;
    string binNum = "";
    // vector<string> subKmers;
};
ostream& operator<<(ostream& out, const data_t& obj);
// map from kmer to vector of indecies where the smer appears in them (used in an Smap)
using Kmap_t = unordered_map<string, vector<int>>;

#endif /* classData_hpp */
