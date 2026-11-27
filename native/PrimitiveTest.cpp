#include "MNA.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
using namespace optical;
int main(){
 MNA m(3);
 auto stamp=[&](){m.source(0,-1,2,1);m.conductance(0,1,.001);m.conductance(1,-1,.001);};
 stamp();if(!m.matrix.prepare())return 1;
 m.matrix.clear();stamp();if(!m.matrix.solve())return 2;
 for(int i=0;i<3;++i)m.v[i]+=m.matrix.x[i];
 if(std::abs(m.v[0]-1)>1e-9||std::abs(m.v[1]-.5)>1e-9)return 3;
 MNA d(3);d.v={1,.5,0};
 auto diode=[&](){d.source(0,-1,2,1);d.conductance(0,1,.001);d.diode(1,-1,1e-12,1);};
 diode();if(!d.matrix.prepare())return 4;
 for(int i=0;i<40;++i){d.matrix.clear();diode();if(!d.matrix.solve())return 5;for(int k=0;k<3;++k)d.v[k]+=std::clamp(d.matrix.x[k],-.1,.1);}
 double res=(1-d.v[1])*.001-1e-12*std::expm1(d.v[1]/thermal)-1e-12*d.v[1];
 if(std::abs(res)>1e-9)return 6;
 std::cout<<"PASS resistor divider "<<m.v[1]<<" V; diode "<<d.v[1]<<" V residual "<<res<<" A\n";
}
