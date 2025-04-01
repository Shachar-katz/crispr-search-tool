//
//  classData.cpp
//  repeatesSearchProjectAltK
//
//  Created by Sarah Katz on 12/29/24.
//

#include <stdio.h>
#include "classData.hpp"

// overloaded stream operator to be able to stream data_t objects
ostream& operator<<(ostream& out, const data_t& obj){
    if (obj.countInFile == 0 && obj.binNum != ""){
        return out << 
         obj.numLines <<
         '\t' << obj.binNum <<
         '\t' << obj.palindromicScore <<
         '\t' << obj.kLen;
    }
    else{
        return out << obj.countInFile << '\t' << obj.numLines;
    }
}
// overloaded for the binsData


