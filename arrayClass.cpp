#include "arrayClass.h"
ostream& operator<<(ostream& out, const Array& obj){
    return out << 
    obj.getRepeatId() <<
    '\t' << obj.getNumSpacers() <<
    '\t' << obj.getArrayLen();
}