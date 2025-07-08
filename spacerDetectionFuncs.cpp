// #include "functions.hpp"
// #include "classData.hpp"
// #include "fileReadClass.hpp"
// #include "arrayClass.h"
// #include "spacerClass.h"

// void scrapeSpacerData(unordered_map<string,Array> globalArrayMap, unordered_map<string,Spacer>& globalSpacerMap){
//     for (auto& [arrayId, array] : globalArrayMap){
//         vector<string> arrayVect = array.getArrayVect();
//         for(int i = 0; i < arrayVect.size(); i++){
//             if (i % 2 == 0){
//                 string currSpacer = arrayVect.at(i);
//                 auto& spacerData = globalSpacerMap[currSpacer];
//                 spacerData.abundance++;
//                 spacerData.associatedArrayId.emplace_back(arrayId);
//                 spacerData.associatedRepeatIdToRepeat[array.getRepeatId()] = array.getRepeat();
//             }
//         }
//     }
// }
