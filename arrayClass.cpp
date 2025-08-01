#include "arrayClass.h"
ostream& operator<<(ostream& out, const Array& obj){
    return out << 
    obj.getRepeatId() <<
    '\t' << obj.getNumSpacers() <<
    '\t' << obj.getArrayLen();
}

ostream& operator<<(ostream& out, const RepeatData& obj){
    return out << 
    obj.repeatID <<
    '\t' << obj.repeatIdx <<
    '\t' << obj.numMissmatches <<
    '\t' << obj.repeatLen;
}