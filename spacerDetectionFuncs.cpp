#include "functions.hpp"
#include "classData.hpp"
#include "fileReadClass.hpp"
#include "arrayClass.h"
#include "spacerClass.h"
using namespace std;

void spacerScraping(const unordered_map<string,Array>& globalArrayMap, unordered_map<string,Spacer>& globalSpacerMap){
    int spacerIdNum = 1;
    for (const auto& [arrayId, array] : globalArrayMap){
        vector<string> arrayVect = array.getArrayVect();
        for(int i = 0; i < arrayVect.size(); i++){
            if (i % 2 == 0){
                string currSpacer = arrayVect.at(i);
                auto& spacerData = globalSpacerMap[currSpacer];
                if (spacerData.abundance == 0)
                {
                    spacerData.spacerId = "S_" + to_string(spacerIdNum);
                    spacerIdNum++;
                }
                spacerData.abundance++;
                // spacerData.associatedArrayId.emplace_back(arrayId);
                // spacerData.associatedRepeatIdToRepeat[array.getRepeatId()] = array.getRepeat();
            }
        }
    }
}
