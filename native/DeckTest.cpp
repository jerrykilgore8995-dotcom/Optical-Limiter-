#include "DeckGenerated.hpp"
#include <iostream>
#include <cmath>
int main(){
 using namespace optical::generated;
 if(nodeCount<100||resistors.size()<100||capacitors.size()<50)return 1;
 for(auto x:resistors)if(!(x.value>0)||!std::isfinite(x.value)||x.a==x.b)return 2;
 for(auto x:capacitors)if(x.value<0||!std::isfinite(x.value)||x.a==x.b)return 3;
 std::cout<<"PASS generated passive deck: "<<nodeCount<<" nodes, "<<resistors.size()<<" resistors, "<<capacitors.size()<<" capacitors; "<<unsupportedCount<<" explicitly pending devices\n";
}
