#pragma once
#include <unordered_map>
#include <iostream>
#include <fstream>
#include "classData.hpp"
#include "fileReadClass.hpp"

// global functions:
string reverseComplement(string seq);
string pickKey(string seq);
template <typename T>
void writeUnorderedMapToFile(const unordered_map<string, T>& map, ofstream& outFS) {
    for (const auto& [key, data] : map) {
        outFS << setw(20) << key << '\t' << data << '\n';

    }
}