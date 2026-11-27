#pragma once
#include "DeviceMath.hpp"
#include "Sparse.hpp"
#include <vector>
#include <stdexcept>
namespace optical {
struct MNA {
    Sparse matrix;
    std::vector<double> v;
    explicit MNA(int n):v(n,0.) {matrix.resize(n);}
    double voltage(int i) const {return i<0?0.:v[i];}
    void residual(int a,int b,double i) {
        if(a>=0)matrix.rhs[a]-=i;
        if(b>=0)matrix.rhs[b]+=i;
    }
    void derivative(int a,int b,int p,int n,double g) {
        matrix.add(a,p,g);matrix.add(a,n,-g);
        matrix.add(b,p,-g);matrix.add(b,n,g);
    }
    void conductance(int a,int b,double g) {
        residual(a,b,(voltage(a)-voltage(b))*g);derivative(a,b,a,b,g);
    }
    void capacitor(int a,int b,double c,double previousV,double dt) {
        residual(a,b,c/dt*(voltage(a)-voltage(b)-previousV));
        derivative(a,b,a,b,c/dt);
    }
    // Trapezoidal companion with accepted previous voltage and current.
    void capacitorTrap(int a,int b,double c,double previousV,double previousI,double dt) {
        const double g=2*c/dt;
        residual(a,b,g*(voltage(a)-voltage(b)-previousV)-previousI);
        derivative(a,b,a,b,g);
    }
    void charge(int a,int b,Jet q,double previousQ,double previousI,double dt) {
        residual(a,b,2*(q.value-previousQ)/dt-previousI);
        derivative(a,b,a,b,2*q.p/dt);
    }
    void source(int a,int b,int k,double volts) {
        residual(a,b,voltage(k));derivative(a,b,k,-1,1.);
        residual(k,-1,voltage(a)-voltage(b)-volts);derivative(k,-1,a,b,1.);
    }
    void controlled(int a,int b,int p,int n,double gain) {
        residual(a,b,gain*(voltage(p)-voltage(n)));derivative(a,b,p,n,gain);
    }
    // Voltage controlled voltage source: branch current plus constraint equation.
    void vcvs(int a,int b,int k,int p,int n,double gain) {
        source(a,b,k,gain*(voltage(p)-voltage(n)));
        derivative(k,-1,p,n,-gain);
    }
    void cccs(int a,int b,int k,double gain) {
        residual(a,b,gain*voltage(k));derivative(a,b,k,-1,gain);
    }
    void ccvs(int a,int b,int k,int control,double gain) {
        source(a,b,k,gain*voltage(control));
        derivative(k,-1,control,-1,-gain);
    }
    void diode(int a,int b,double is,double n) {
        const auto j=junction({voltage(a)-voltage(b),1,0},n*thermal,is);
        residual(a,b,j.value+1e-12*(voltage(a)-voltage(b)));
        derivative(a,b,a,b,j.p+1e-12);
    }
    void terminal(int row,int c,int b,int e,int sign,Jet j) {
        if(row>=0)matrix.rhs[row]-=sign*j.value;
        matrix.add(row,c,-j.q);matrix.add(row,b,j.p+j.q);matrix.add(row,e,-j.p);
    }
    void jfet(int d,int g,int s,int sign,double vto,double beta,double is,double lambda=0) {
        const Jet gs{sign*(voltage(g)-voltage(s)),1,0},gd{sign*(voltage(g)-voltage(d)),0,1};
        const Jet ds=gs-gd;const bool reverse=ds.value<0;
        const Jet u=(reverse?gd:gs)-Jet{vto};
        const Jet vd=ds*(reverse?-1:1);
        const Jet b=Jet{beta}+vd*(beta*lambda);
        const Jet ids=u.value<=0?Jet{0}:b*(vd.value>=u.value?u*u:vd*(u*2-vd));
        const Jet id=ids*(reverse?-1:1);
        const Jet igs=junction(gs,thermal,is)+gs*1e-12,igd=junction(gd,thermal,is)+gd*1e-12;
        terminal(d,d,g,s,sign,id-igd);terminal(g,d,g,s,sign,igs+igd);terminal(s,d,g,s,sign,(id+igs)*-1);
    }
    // DC Gummel-Poon transport. Dynamic charge and series resistances follow in deck integration.
    struct BJT {double is,bf,nf,br,nr,vaf,var,ikf,ikr,ise,ne,isc,nc;};
    void bjt(int c,int b,int e,int sign,const BJT& m) {
        const Jet be{sign*(voltage(b)-voltage(e)),1,0},bc{sign*(voltage(b)-voltage(c)),0,1};
        const Jet f=junction(be,m.nf*thermal,m.is),r=junction(bc,m.nr*thermal,m.is);
        Jet base=Jet{1}-be*(m.var?1/m.var:0)-bc*(m.vaf?1/m.vaf:0);
        if(m.ikf||m.ikr) {
            Jet k=Jet{1}+(f*(m.ikf?1/m.ikf:0)+r*(m.ikr?1/m.ikr:0))*4;
            if(k.value<1e-12)k={1e-12};
            base=base/(Jet{.5}+squareRoot(k)*.5);
        }
        const Jet t=(f-r)*base,ibf=f*(1/m.bf)+junction(be,m.ne*thermal,m.ise)+be*1e-12,ibr=r*(1/m.br)+junction(bc,m.nc*thermal,m.isc)+bc*1e-12;
        terminal(c,c,b,e,sign,t-ibr);terminal(b,c,b,e,sign,ibf+ibr);terminal(e,c,b,e,sign,(t+ibf)*-1);
    }
};
}
