#pragma once
#include <iostream>
#include <string>
using namespace std;

class Spacer{
public:
    string spacerId;
    string spacer;
    unordered_map<string,string> associatedRepeatIdToRepeat;
    vector<string> associatedArrayId;
    int abundance = 0;
};

// // // ostream& operator<<(ostream& out, const Array& obj);

// // // ostream& operator<<(ostream& out, const RepeatData& obj);