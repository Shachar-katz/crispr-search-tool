#include <iostream>
#include <string>
using namespace std;

class Array{
private:
    int startIndex;
    int endIndex;
    string repeat;
    string array;
    int arrayLen;
    bool isOpen;
public:
    Array(int idxInLine, string repeat)
    {
        startIndex = idxInLine;
        this->repeat = repeat;
        isOpen = true;
    }
    ~Array(){}
    bool isArrayOpen(){ return this->isOpen; }
    int getStartIndex() { return this->startIndex; }
    int getEndIndex() { return this->endIndex; }
    void closeArray(int endIdxInLine, string line) {
        endIndex = endIdxInLine;
        isOpen = false;
        arrayLen = endIndex - startIndex;
        array = line.substr(startIndex, arrayLen);
    }

};