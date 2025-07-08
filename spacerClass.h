#pragma once
#include <iostream>
#include <string>
using namespace std;

class Spacer{
public:
    string spacerId;
    int abundance = 0;
};

ostream& operator<<(ostream& out, const Spacer& obj);

// // // ostream& operator<<(ostream& out, const RepeatData& obj);