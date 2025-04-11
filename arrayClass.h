#include <iostream>
#include <string>
using namespace std;

class Array{
private:
    string repeat;
    vector<int> inLineCoordinatesVect;
    vector<int> inArrayRepeatCoordinatesVect;
    vector<string> array;
    int arrayLen;
    bool isOpen; // might not need 
public:
    Array(){}
    Array(string repeat)
    {
        this->repeat = repeat;
        isOpen = true; // might not need 
    }
    ~Array(){}
    void openArray(string repeat, int startIdxInLine)
    {
        this->repeat = repeat;
        isOpen = true; // might not need 
        addRepeat(startIdxInLine);
    }
    void addRepeat(int startIdxInLine)
    {   
        int endIdxInLine = startIdxInLine + repeat.length();
        inLineCoordinatesVect.push_back(startIdxInLine);
        inLineCoordinatesVect.push_back(endIdxInLine);
    }
    bool stillValid(int startIdxInLine, int maxValidSpacer){
        int spacerLen = startIdxInLine - inLineCoordinatesVect.back();
        if (spacerLen <= maxValidSpacer){ return true; }
        return false;
    }
    bool isArrayOpen(){ return this->isOpen; } // might not need 
    string getRepeat() { return repeat; }
    bool closeArray(string line) 
    {   
        isOpen = false; // might not need 
        if (inLineCoordinatesVect.size() < 4) { return false; }
        arrayLen = inLineCoordinatesVect.back() - inLineCoordinatesVect.at(0);
        for (int i = 0; i < inLineCoordinatesVect.size() - 1; i++)
        {
            int segmentLen = inLineCoordinatesVect.at(i + 1) - inLineCoordinatesVect.at(i);
            string arraySegment = line.substr(inLineCoordinatesVect.at(i), segmentLen);
            array.push_back(arraySegment);
            // potentially populate repeat in array coordinate vect
        }
        return true;
    }

};

class LineArrayHandler{
private:
    string line;
    int maxAllowedSpacer;
    vector<Array> lineArrayVect;
    bool activeArray = false;
    Array tempArray;
public:
    LineArrayHandler(const string& line, int maxAllowedSpacer)
    {
        this->line = line;
        this->maxAllowedSpacer = maxAllowedSpacer;
    }
    ~LineArrayHandler(){}
    void manageState(string repeat, int startIdxInLine)
    {
        if(!activeArray){
            tempArray.openArray(repeat, startIdxInLine);
            activeArray = true;
        }
        else if(tempArray.getRepeat() == repeat){
            this->expandArray(startIdxInLine);
        }
        else if(tempArray.getRepeat() != repeat){
            this->uploadArray(); // should I use return value for anything?
            tempArray.openArray(repeat, startIdxInLine);
            activeArray = true;
        }
    }
    bool expandArray(int startIdxInLine)
    {
        if (tempArray.stillValid(startIdxInLine, maxAllowedSpacer)){
            tempArray.addRepeat(startIdxInLine);
            return true;
        }
        return uploadArray();
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