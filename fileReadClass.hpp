//
//  fileReadClass.hpp
//  repeatesSearchProjectAltK
//
//  Created by Sarah Katz on 2/7/25.
//

#ifndef fileReadClass_h
#define fileReadClass_h

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <memory>

using namespace std;

enum class FileType {
    TXT,
    FASTA,
    FASTQ_SINGLE,
};

// helper parser
FileType parseFileType(const string &fileTypeStr);

// abstract base class for custome get line 
class customGetline {
public:
    virtual ~customGetline() = default;
    virtual bool getNextRead(string &outLine) = 0;
};

// getline txt class: validates file and returns a bool getline
class txtGetline : public customGetline {
private:
    ifstream &inFile;
public:
    explicit txtGetline(ifstream &fileStream) : inFile(fileStream) {}
    bool getNextRead(string &outLine) override {
        outLine.clear();
        if (!inFile.good()) {
            return false;
        }
        return static_cast<bool>(getline(inFile, outLine));
    }
};

// fatsa getline: takes only lines after > line and concatenates with all the reads until the next >
class fastaGetline : public customGetline {
private:
    ifstream &inFile;
public:
    explicit fastaGetline(ifstream &fileStream) : inFile(fileStream) {}
    bool getNextRead(string &outLine) override {
        outLine.clear();
        if (!inFile.good()) {
            return false;
        } // test if file is good
        string tempLine;
        while (getline(inFile, tempLine)) {
            if (!tempLine.empty() && tempLine[0] == '>') {
                break;
            }
        } // get to the next header
        if (!inFile.good() || tempLine.empty() || tempLine[0] != '>') {
            return false;
        } // if file is over/ there is no more header return false

        outLine.clear(); //clear the line we are updating
        streampos savedPosition = inFile.tellg(); //update position of pointer in file

        while (getline(inFile, tempLine)) {
            if (!tempLine.empty() && tempLine[0] == '>') {
                // Found the header of the *next* record
                // Rewind so next call can start there
                inFile.seekg(savedPosition);
                break;
            }
            outLine += tempLine; // add up the reads in the segment
            savedPosition = inFile.tellg();
        }
        return !outLine.empty();
    }
};

// fastq single getline: takes every second line out of 4 
class fastqSingleGetline : public customGetline {
private:
    ifstream &inFile;
public:
    explicit fastqSingleGetline(ifstream &fileStream) : inFile(fileStream) {}
    bool getNextRead(string &outLine) override {
        outLine.clear();
        if (!inFile.good()) {
            return false;
        }
        string discard;
        string seqLine;

        if (!getline(inFile, discard)) return false;
        if (!getline(inFile, seqLine)) return false;
        if (!getline(inFile, discard)) return false;
        if (!getline(inFile, discard)) return false;

        outLine = seqLine;
        return true;
    }
};

// helper function to generate pointer to right kind of reader (unique pointer self memory manages)
unique_ptr<customGetline> createReader(FileType fileType, ifstream &inFile);

// manages file validation, file opening, and generating the right reader, then creates a getline for that reader
class MultiFormatFileReader {
private:
    ifstream inFile;
    unique_ptr<customGetline> reader;
    FileType fileType;
    
public:
    // the manager opens file, parses the file type, and creates the correct reader
    MultiFormatFileReader(const string &filePath, const string &fileTypeStr) {
        inFile.open(filePath);
        if (!inFile.is_open()) {
            throw runtime_error("Could not open file: " + filePath);
        }
        fileType = parseFileType(fileTypeStr);
        reader = createReader(fileType, inFile);
        if (!reader) {
            throw runtime_error("Failed to create reader for file type, most likely file type arg was entered wrong.");
        }
    }
    bool getNextLine(string &line) {
            return reader->getNextRead(line);
    }
    FileType getFileType(){
        return this->fileType;
    }
    
};
#endif /* fileReadClass_h */
