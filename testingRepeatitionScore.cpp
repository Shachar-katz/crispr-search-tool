// #include "functions.hpp"
// #include "classData.hpp"
// #include "fileReadClass.hpp"




// // this function populates the smap
// void findSeedPattern(string line, unordered_map<string,vector<int>>& singleLineMapSeedKToIdx, int seedK){
//     for (int j = 0; j <= (line.size() - seedK); j++){
//         string key = line.substr(j,seedK);
//         if (key.size() < seedK){
//             break; // throw error
//         }
//         singleLineMapSeedKToIdx[key].emplace_back(j);
//     }
// }

// // this function recives start and end location of 2 Kmers in expansion and 
// // varifies that they are not overlapping
// bool notOverlapping(int startIdx, int idxStartCompare, int idxEnd, int idxEndCompare, int minLegitimateSpacer)
// {
//     int spacing;
//     if (idxEnd < idxStartCompare) {
//         spacing = idxStartCompare - idxEnd -1;
        
//     } 
//     else if (idxEndCompare < startIdx) {
//         spacing = startIdx - idxEndCompare - 1;
//     }
//     else {
//         return false;
//     }
//     return (spacing >= minLegitimateSpacer);
// }

// // this function populates the unique kmer set
// int expandSeedTesting(const string& line, 
//                       string smer, 
//                       int startIdx, 
//                       vector<int> smerIdxVect, 
//                       int minK, 
//                       unordered_map<int, string>& uniqueKmersInLine,
//                       int minLegitimateSpacer, 
//                       int maxK, 
//                       int horizion){
//     // Track unique K-mers for this particular smer on this particular line
//     // create range
//     int lowerBound = startIdx - horizion;
//     int upperBound = startIdx + horizion;

//     auto lowerIt = lower_bound(smerIdxVect.begin(), smerIdxVect.end(), lowerBound);
//     auto upperIt = upper_bound(smerIdxVect.begin(), smerIdxVect.end(), upperBound);

//     string bestK = "";
//     int bestKPos = -1;

//     // cout << "expanded for s = " << smer << " at location i = " << startIdx;

//     // then we iterate over all the other indecies comparing the smer appearance there to our "Potential kmer"
//     for(auto it = lowerIt; it != upperIt; it++){
//         if (*it == startIdx){ continue; }
//         // cout << " Against " << *it << endl;

//         // we set the indecies of the start and end at our potential kmer and our comparision
//         int startIdxCopy = startIdx;
//         int idxEnd = startIdx + smer.length() - 1;
        
//         int idxStartCompare = *it;
//         int idxEndCompare = idxStartCompare + smer.length() - 1;
        
//         // we set flags = true to represent are the 2 occurances equal at the start and end.
//         bool equalAtStart = true;
//         bool equalAtEnd = true;
//         // the current Kmer is a size of an smer
//         int kmerLen = smer.length();
//         // while the 2 locations are equal at the start or end and the indecies are not overlapping:
//         while((equalAtStart || equalAtEnd) && 
//                 notOverlapping(startIdxCopy, idxStartCompare, idxEnd, idxEndCompare, minLegitimateSpacer) &&
//                 kmerLen < maxK){
//             // as long as the indecies would valid at the start if we decremented and they are still equal at the start:
//             if (equalAtStart && startIdxCopy > 0 && idxStartCompare > 0){
//                 // if the "next" nucleotid is equal on both occurances we decrement the position and continue expending
//                 if (line[startIdxCopy - 1] == line[idxStartCompare - 1]){
//                     startIdxCopy-- ;
//                     idxStartCompare-- ;
//                 }
//                 // else we turn the equal at start flag to false.
//                 else{
//                     equalAtStart = false;
//                 }
//             }
//             else{
//                 equalAtStart = false;
//             }
//             // after one end expansion if they overlapp we want to stop
//             if(!notOverlapping(startIdxCopy,idxStartCompare,idxEnd,idxEndCompare,minLegitimateSpacer)) {break;}
//             // as long as the indecies would be valid at the end if we incremented and they are still equal at the end:
//             if (equalAtEnd && (idxEnd + 1) < line.length() && (idxEndCompare + 1) < line.length()){
//                 // if the "next" nucleotid is equal on both occurances we increment the position and continue expending
//                 if (line[idxEnd + 1] == line[idxEndCompare + 1]){
//                     idxEnd++ ;
//                     idxEndCompare++ ;
//                 }
//                 // else we turn the equal at end flag to false.
//                 else{
//                     equalAtEnd = false;
//                 }
//             }
//             else{
//                 equalAtEnd = false;
//             }
//             kmerLen = idxEnd - startIdxCopy + 1;    
//         } 
//         // if they reached a point where either Kmers overlap, or hit any of the bounds, throw them out
//         // if we are strict we also throw out the entire line
//         if (!notOverlapping(startIdxCopy,idxStartCompare,idxEnd,idxEndCompare,minLegitimateSpacer) || 
//             idxEnd >= (line.length() - 1) || 
//             idxEndCompare >= (line.length() - 1) || 
//             startIdxCopy <= 0 || 
//             idxStartCompare <= 0 ||
//             kmerLen >= maxK){
//                 continue; 
//         }
//         // now if we reached this part of the code we are sure we have reached max expansion:
//         // 1) the nucleotids dont match anymore
//         // 2) they haven't reached a point of hitting boundaries or overlapping  (overlap defined as include spacer)
        
//         // if this potential kmer follows requirements :
//         // (is indeed a kmer and is longer then the best kmer so far)
//         // we make it the best kmer
//         if(kmerLen >= minK && kmerLen > bestK.length()){
//             string kmer = line.substr(startIdxCopy, kmerLen);
//             // cout << "bestK for s = " << smer << " at location i = " << startIdxCopy << " resulted in k= " << kmer << endl; //db
//             bestK = kmer;
//             bestKPos = startIdxCopy;
//         }   
//     }
//     // once we have compared every smer to all the other Smers:
//     // We add the best match for each instance (a kmer) to the unique Kmers in line Set.
//     if(bestKPos != -1){
//         auto& repeat = uniqueKmersInLine[bestKPos];
//         repeat = bestK;
//         return 1;
//     }
//     return 0;
// }

// void generateRepeatition(const string& line,
//                          int segmentSize,
//                          int seedK, 
//                          int minK,
//                          int maxK,
//                          int minLegitimateSpacer,
//                          int horizon,
//                          int smoothingWindow,
//                          const unordered_map<string,vector<int>>& singleLineMapSeedKToIdx, 
//                          unordered_map<int,double>& inLineSmoothRepetition)
// {
//     unordered_map<int,string> posToKmerInLine;
//     vector<int> inLineSegments;
//     vector<double> inLineRepetitionScores;
//     unordered_map<int,double> inLineSegmentToRepetition;

//     int currSegment = 0;
//     double currRepitionScore = 0;

//     for (int i = 0; i <= (line.length() - seedK) ; i++){
//         string smer = line.substr(i,seedK);
//         auto& idxs = singleLineMapSeedKToIdx.at(smer);
//         if (idxs.size() > 1){
//             currRepitionScore += expandSeedTesting(line, smer, i, idxs, minK, posToKmerInLine, minLegitimateSpacer, maxK, horizon);
//         }
//         if (i % segmentSize == 0){
//             inLineSegments.emplace_back(currSegment);
//             inLineRepetitionScores.emplace_back(currRepitionScore / segmentSize);
//             currSegment = i; 
//             currRepitionScore = 0.0;
//         }
//     }

//     if (inLineSegments.size() != inLineRepetitionScores.size()) { cout << "error" << endl;/* throw error*/}
//     const int numSegments = static_cast<int>(inLineRepetitionScores.size());

//     for (int j = 0; j < numSegments; j++)
//     {
//         int lowIdx = std::max(j - smoothingWindow, 0);
//         int highIdx = std::min(j + smoothingWindow, numSegments - 1);
//         double smoothedScore = 0.0;
//         double sumA = 0.0;
//         for (int i = lowIdx; i <= highIdx; i++){
//             int aI = 2 * smoothingWindow - abs(i - j);
//             double weightI = static_cast<double>(aI);
//             smoothedScore += weightI * inLineRepetitionScores[i];
//             sumA += aI;
//         }
//         int pos = inLineSegments.at(j);
//         inLineSmoothRepetition[pos] = smoothedScore / sumA;
//     }
// }

// void generateRepeatitionOverLines(MultiFormatFileReader& fileReader, 
//                      int seedK, 
//                      int minK, 
//                      int minLegitimateSpacer,
//                      int maxLegitimateSpacer,
//                      int horizon,
//                      ofstream& repitionChart,
//                      int interval,
//                      int maxK)
// {
//     // line variable temporerally holds the reads
//     string line;
//     // statistics Vars
//     int progressCounter = 0;
//     int numReadsWithRepeats = 0;
//     int faultyLine = 0;
//     // loop over every read
//     while (fileReader.getNextLine(line)) {
//         // statistics and progress managment:
//         progressCounter++;
//         if (progressCounter % interval == 0){
//             cout << "Procession line: " << progressCounter << endl;
//         }
        
//         // every line we create and populate an smap that maps from an smer to vect of indecies in the line.
//         unordered_map<string,vector<int>> singleLineMapSeedKToIdx;
//         findSeedPattern(line, singleLineMapSeedKToIdx, seedK);
//         int segmentSize = 100;
//         int smoothingWindow = 2;
//         // we also create an empty map of unique repeats -> set of their positions in this line.
//         unordered_map<int,double> inLineSmoothRepetition;
//         generateRepeatition(line, segmentSize, seedK, minK, maxK, minLegitimateSpacer, horizon, smoothingWindow, singleLineMapSeedKToIdx, inLineSmoothRepetition);

//         for (const auto& [location, repetiton] : inLineSmoothRepetition) {
//             repitionChart << '\t' << progressCounter << '\t' << location << '\t' << repetiton << endl;;
//         }
//     }
// }


// int main(){
//     // open log file
//     ofstream repitionChart;
//     string repitionChartFile = "/Users/sarahkatz/relman_lab/L1_V5/repeatition_chart_bad_lines_test";
//     repitionChart.open(repitionChartFile);
//     if (!repitionChart.is_open()){
//          cerr << "Error: Could not open log output file." << endl;
//          return -1;
//     }
//     repitionChart << "read_num" << '\t' << "block_coordinations" << '\t' << "repeatition" << endl;

//     string inputFile = "/Users/sarahkatz/relman_lab/PacBio/new_bad_lines.txt";
//     MultiFormatFileReader fileReader(inputFile, "txt");

//      // calculate seedK and horizon
//     int seedK = 10;
//     int minK = 20;
//     int maxK = 70;
//     int minLegitSpace = 10;
//     int maxLegitSpace = 100;
//     int interval = 1;
//     int horizon = (maxK + maxLegitSpace) * 4 + maxK;
//     int L = 150;
//      // populate the global Kmer map (identify repeats)
//     generateRepeatitionOverLines(fileReader, seedK, minK, minLegitSpace, maxLegitSpace, horizon, repitionChart, interval, maxK);    
//     repitionChart.close();
//     return 0;
// }
// // int main(){
// //     // open log file
// //     ofstream repitionChart;
// //     string repitionChartFile = "/Users/sarahkatz/relman_lab/L1_V4/repeatition_chart_bad_lines";
// //     repitionChart.open(repitionChartFile);
// //     if (!repitionChart.is_open()){
// //          cerr << "Error: Could not open log output file." << endl;
// //          return -1;
// //     }
// //     repitionChart << "read_num" << '\t' << "block_coordinations" << '\t' << "repeatition" << endl;

// //     string inputFile = "/Users/sarahkatz/relman_lab/PacBio/badLines.fa";
// //     MultiFormatFileReader fileReader(inputFile, "fasta");

// //      // calculate seedK and horizon
// //     int seedK = 10;
// //     int minK = 20;
// //     int maxK = 70;
// //     int minLegitSpace = 10;
// //     int maxLegitSpace = 100;
// //     int interval = 1;
// //     int horizon = (maxK + maxLegitSpace) * 4 + maxK;
// //     int L = 150;
// //      // populate the global Kmer map (identify repeats)
// //     generateRepeatition(fileReader, seedK, minK, minLegitSpace, maxLegitSpace, horizon, repitionChart, interval, maxK);    
// //     repitionChart.close();
// //     return 0;
// // }
