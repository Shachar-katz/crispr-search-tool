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
        setw(10) << obj.numLines <<
        setw(10) << '\t' << obj.binNum <<
        setw(10) << '\t' << obj.palindromicScore <<
        setw(10) << '\t' << obj.KLen << endl;
    }
    else{
        return out << setw(30) << obj.countInFile << '\t' << obj.numLines << '\n';
    }
}
// overloaded for the binsData


