#pragma once
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include "classData.hpp"
#include "fileReadClass.hpp"

// global functions:
string reverseComplement(string seq);
string pickKey(string seq);
template <typename T>
void writeUnorderedMapToFile(const unordered_map<string, T>& map, ofstream& outFS) {
    for (const auto& [key, data] : map) {
        outFS << key << '\t' << data << endl;

    }
}
struct File {
    string identifier;
    string filePath;
};

// New: function to read identifier table
vector<File> readIdentifierTable(const string& baseDirectory, const string& tableFile);