#include "fileReadClass.hpp"

// parases file type
FileType parseFileType(const string &fileTypeStr) {
    if (fileTypeStr == "txt") {
        return FileType::TXT;
    } else if (fileTypeStr == "fasta") {
        return FileType::FASTA;
    } else if (fileTypeStr == "fastq_single" || fileTypeStr == "fastq_dual") {
        return FileType::FASTQ_SINGLE;
    } else {
        throw runtime_error("Unsupported file type: " + fileTypeStr);
    }
}

// helper function to generate pointer to right kind of reader (unique pointer self memory manages)
unique_ptr<customGetline> createReader(FileType fileType, ifstream &inFile) {
    switch (fileType) {
        case FileType::TXT:
            return make_unique<txtGetline>(inFile);
        case FileType::FASTA:
            return make_unique<fastaGetline>(inFile);
        case FileType::FASTQ_SINGLE:
            return make_unique<fastqSingleGetline>(inFile);
    }
    throw runtime_error("Unknown FileType enum encountered.");
}