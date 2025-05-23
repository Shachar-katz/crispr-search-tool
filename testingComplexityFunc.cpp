#include "functions.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"

void generateComplexities(MultiFormatFileReader& fileReader, 
                     int seedK, 
                     int segmentsSize, 
                     ofstream& complexityChart)
    {
    // line variable temporerally holds the reads
    string line;
    // statistics Vars
    int progressCounter = 0;
    // loop over every read
    while (fileReader.getNextLine(line)) {
        // statistics and progress managment:
        progressCounter++;
        cout << "Procession line: " << progressCounter << endl;
        
        // every line we create and populate an smap that maps from an smer to vect of indecies in the line.
        unordered_map<int,double> inLineSegmentToComplexity;

        for (int i = 0; i <= (line.length() - segmentsSize) ; i += segmentsSize){
            string segment = line.substr(i,segmentsSize);
            unordered_set<string> uniqueSmersInSegment;
            findSmerSet(segment, uniqueSmersInSegment, seedK);
            double complexity = static_cast<double>(uniqueSmersInSegment.size()) / static_cast<double>(segmentsSize);
            inLineSegmentToComplexity[i] = complexity;
        }

        for (const auto& [location, complexity] : inLineSegmentToComplexity) {
            complexityChart << '\t' << progressCounter << '\t' << location << '\t' << complexity << endl;;
        }
    }
}

int main(){
    // open log file
    ofstream complexityChart;
    string complexityChartFile = "/Users/sarahkatz/relman_lab/L1_V4/complexity_chart_array_lines";
    complexityChart.open(complexityChartFile);
    if (!complexityChart.is_open()){
         cerr << "Error: Could not open log output file." << endl;
         return -1;
    }
    complexityChart << "read_num" << '\t' << "block_coordinations" << '\t' << "complexity" << endl;

    string inputFile = "/Users/sarahkatz/relman_lab/PacBio/ideating.txt";
    MultiFormatFileReader fileReader(inputFile, "txt");

     // calculate seedK and horizon
    int seedK = 10;
    int L = 150;
     // populate the global Kmer map (identify repeats)
    generateComplexities(fileReader, seedK, L, complexityChart);    
    complexityChart.close();
    return 0;
}
