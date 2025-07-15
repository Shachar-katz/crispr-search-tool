#pragma once
#include <iostream>
#include <string>
#include <vector>
using namespace std;

class Array{
private:
    string repeat;
    string repeatId;
    vector<int> inLineCoordinatesVect;
    vector<int> inArrayRepeatCoordinatesVect; // might not need?
    vector<string> array;
    vector <int> repeatIdxToNumMissmatches;
    // add spacer object?
    int arrayLen = 0; // might not need?
    int numSpacers;
public:
    Array(){}
    Array(string repeat)
    {
        this->repeat = repeat;
    }
    ~Array(){}
    void openArray(string repeat, int startIdxInLine, string repeatId, int numMissmatches)
    {
        this->repeat = repeat;
        this->repeatId = repeatId;
        addRepeat(startIdxInLine, numMissmatches);
    }
    void addRepeat(int startIdxInLine, int numMissmatches)
    {   
        int endIdxInLine = startIdxInLine + repeat.length();
        inLineCoordinatesVect.push_back(startIdxInLine);
        inLineCoordinatesVect.push_back(endIdxInLine);
        repeatIdxToNumMissmatches.push_back(numMissmatches);
    }
    bool stillValid(int startIdxInLine, int maxValidSpacer, int minValidSpacer){
        int spacerLen = startIdxInLine - inLineCoordinatesVect.back();
        if (spacerLen <= maxValidSpacer && startIdxInLine > (inLineCoordinatesVect.back() + minValidSpacer)){ return true; }
        return false;
    }
    string getRepeat() const { return repeat; }
    string getRepeatId() const { return repeatId; }
    int getNumSpacers() const { return numSpacers; }
    int getArrayLen() const { return arrayLen; }
    vector<string> getArrayVect() const { return this->array; }
    string getArrayStr() 
    {
        string array;
        for (int i = 0; i < this->array.size(); i++)
        {
            array += this->array[i] + "|";
        }
        return array;
    }
    vector<int> getErrorVect() const { return repeatIdxToNumMissmatches; }
    bool closeArray(string line){
        if (inLineCoordinatesVect.size() < 4) { return false; }
        arrayLen = inLineCoordinatesVect.back() - inLineCoordinatesVect[0];
        for (int i = 0; i < inLineCoordinatesVect.size() - 1; i++)
        {
            int segmentLen = inLineCoordinatesVect[i + 1] - inLineCoordinatesVect[i];
            string arraySegment = line.substr(inLineCoordinatesVect.at(i), segmentLen);
            array.push_back(arraySegment);
            // potentially populate repeat in array coordinate vect
        }
        // calculate num spacers
        int numRepeats = inLineCoordinatesVect.size() / 2;
        numSpacers = numRepeats - 1;
        return true;
    }


};

class LineArrayHandler{
private:
    string line;
    int maxAllowedSpacer;
    int minAllowedSpacer;
    vector<Array> lineArrayVect;
    bool activeArray = false;
    Array tempArray;
public:
    LineArrayHandler(const string& line, int maxAllowedSpacer, int minAllowedSpacer)
    {
        this->line = line;
        this->maxAllowedSpacer = maxAllowedSpacer;
        this->minAllowedSpacer = minAllowedSpacer;
    }
    ~LineArrayHandler(){}
    void manageState(string repeat, int startIdxInLine, string repeatId, int numMissmatches)
    {
        if(!activeArray){
            tempArray.openArray(repeat, startIdxInLine, repeatId, numMissmatches);
            activeArray = true;
        }
        else if(tempArray.getRepeat() == repeat){
            bool wasExpanded = this->expandArray(startIdxInLine, numMissmatches);
            if (!wasExpanded){
                tempArray.openArray(repeat, startIdxInLine, repeatId, numMissmatches);
                activeArray = true;
            }
        }
        else if(tempArray.getRepeat() != repeat){
            this->uploadArray(); // should I use return value for anything?
            tempArray.openArray(repeat, startIdxInLine, repeatId, numMissmatches);
            activeArray = true;
        }
    }
    bool expandArray(int startIdxInLine, int numMissmatches)
    {
        if (tempArray.stillValid(startIdxInLine, maxAllowedSpacer, minAllowedSpacer)){
            tempArray.addRepeat(startIdxInLine, numMissmatches);
            return true;
        }
        bool uploaded = uploadArray();
        return false;
    }
    bool uploadArray()
    {
        bool validArray = tempArray.closeArray(line);
        if (validArray){ 
            lineArrayVect.push_back(tempArray);
        }
        tempArray = Array();
        activeArray = false;
        return validArray;
    }
    bool isActive() { return activeArray; }
    vector<Array> getLineArrayVect(){ return lineArrayVect; }
    bool noArrays() {
        if (lineArrayVect.empty()){
            return true;
        }
        return false;
    }
    
};

ostream& operator<<(ostream& out, const Array& obj);

class RepeatData{
public:
    string arrayID;
    string repeatID;
    int repeatIdx;
    int numMissmatches;
    int repeatLen;
};

ostream& operator<<(ostream& out, const RepeatData& obj);