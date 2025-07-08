//
//  functions.cpp
//
//  Created by Shachar Katz on 10/23/24.
//

// #include "functions.hpp"
#include "globalFunctions.h"


// global functions bellow:

bool isInputFileValid(ifstream& inFS, string fileName) {
    if (!inFS.is_open()) {
        cerr << "Error: Could not open the file: " <<  fileName << endl << "It might not exist or be corrupted." << endl;
        return false;
    }
    // Check if the file is empty
    inFS.seekg(0, ios::end);
    if (inFS.tellg() == 0) {
        cerr << "Error: The following file is empty: " << fileName << endl;
        return false;
    }
    inFS.clear();
    inFS.seekg(0, ios::beg);
    return true;
}

string reverseComplement(string seq) {
    string result = seq;
    int N = seq.length();
    for (int i = 0; i < N; i++) {
        char nucleaotid = seq[N - i - 1];
        char reverse;
        switch(nucleaotid)
        {
            case 'A': reverse = 'T'; break;
            case 'G': reverse = 'C'; break;
            case 'C': reverse = 'G'; break;
            case 'T': reverse = 'A'; break;
            default:  reverse = nucleaotid;
        }
        result[i] = reverse;
    }
  return(result);
}

string pickKey(string seq) {
    string reverse = reverseComplement(seq);
    return (reverse > seq) ? reverse : seq;
}

vector<string> readFileList(const string& listFile) {
    vector<string> files;
    ifstream inFile(listFile);
    string line;
    
    if (!inFile.is_open()) {
        cerr << "Error: Could not open file list: " << listFile << endl;
        return files;
    }
    
    while (getline(inFile, line)) {
        // Remove whitespace and skip empty lines
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (!line.empty()) {
            files.push_back(line);
        }
    }
    
    inFile.close();
    return files;
}