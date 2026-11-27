#pragma once
#include "DeviceMath.hpp"
#include "Sparse.hpp"
#include <vector>
#include <stdexcept>
namespace optical {
struct MNA {
    Sparse matrix;
    std::vector<double> v;
    explicit MNA(int n, int fixed=0):v(n+fixed,0.), unknowns(n) {matrix.resize(n); companions.reserve(128);}
    double voltage(int i) const {return i<0?0.:v[i];}
    void residual(int a,int b,double i) {
        if(a>=0 && a<unknowns)matrix.rhs[a]-=i;
        if(b>=0 && b<unknowns)matrix.rhs[b]+=i;
    }
    void derivative(int a,int b,int p,int n,double g, bool variable=false) {
        const int plus=p<unknowns?p:-1, minus=n<unknowns?n:-1;
        if (a<unknowns) {matrix.add(a,plus,g,variable);matrix.add(a,minus,-g,variable);}
        if (b<unknowns) {matrix.add(b,plus,-g,variable);matrix.add(b,minus,g,variable);}
    }
    void conductance(int a,int b,double g,bool variable=false) {
        residual(a,b,(voltage(a)-voltage(b))*g);derivative(a,b,a,b,g,variable);
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
    template<class F> void behavior(int a,int b,int k,std::initializer_list<int> refs,F fn) {
        auto w=[&](int n){return voltage(n);};
        const double value=fn(w);
        if(k>=0)source(a,b,k,value);else residual(a,b,value);
        for(int n:refs) {
            const double h=1e-7*std::max(1.,std::abs(voltage(n)));
            auto hi=[&](int p){return voltage(p)+(p==n?h:0.);};
            auto lo=[&](int p){return voltage(p)-(p==n?h:0.);};
            const double g=(fn(hi)-fn(lo))/(2*h);
            if(k>=0)derivative(k,-1,n,-1,-g,true);else derivative(a,b,n,-1,g,true);
        }
    }
    void diode(int a,int b,double is,double n) {
        const auto j=junction({voltage(a)-voltage(b),1,0},n*thermal,is);
        residual(a,b,j.value+1e-12*(voltage(a)-voltage(b)));
        derivative(a,b,a,b,j.p+1e-12,true);
    }
    void terminal(int row,int c,int b,int e,int sign,Jet j) {
        residual(row,-1,sign*j.value);
        derivative(row,-1,c,-1,-j.q,true);
        derivative(row,-1,b,-1,j.p+j.q,true);
        derivative(row,-1,e,-1,-j.p,true);
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
    // Each capacitor or semiconductor junction keeps its last accepted
    // charge and companion current. Failed Newton trials never advance them.
    struct Companion { double charge = 0, current = 0, trial = 0; };
    void beginCompanions(double coefficient) { inverse = coefficient; nextCompanion = 0; }
    double companion(double q) {
        if (nextCompanion == companions.size()) companions.emplace_back();
        auto& h = companions[nextCompanion++];
        h.trial = q;
        return inverse > 0 ? inverse * (q - h.charge) - h.current : 0;
    }
    void endCompanions() {}
    void acceptCompanions(bool dc) {
        for (auto& h : companions) {
            h.current = dc ? 0 : inverse * (h.trial - h.charge) - h.current;
            h.charge = h.trial;
        }
    }
    void capacitive(int a, int b, double c) {
        residual(a, b, companion(c * (voltage(a) - voltage(b))));
        derivative(a, b, a, b, inverse * c);
    }

    struct Diode {
        double is, n, bv, ibv, cjo, vj, m, fc, tt;
    };
    void diodeDevice(int a, int b, Diode p) {
        const Jet volts{voltage(a) - voltage(b), 1};
        const double vt = p.n * thermal;
        const double breakdown = p.bv - vt * std::log(p.ibv / p.is);
        const auto diffusion = p.bv > 0 && volts.value < -breakdown
            ? junction(Jet{-breakdown} - volts, vt, p.is) * -1 + Jet{-p.is}
            : junction(volts, vt, p.is);
        double i = diffusion.value + 1e-12 * volts.value;
        double g = diffusion.p + 1e-12;
        if (p.cjo != 0 || p.tt != 0) {
            const auto q = depletion(volts, p.cjo, p.vj, p.m, p.fc) + diffusion * p.tt;
            i += companion(q.value);
            g += inverse * q.p;
        }
        residual(a, b, i);
        derivative(a, b, a, b, g, true);
    }

    struct BJT {
        double is, bf, nf, br, nr, vaf, var, ikf, ikr, ise, ne, isc, nc;
        double nk, cje, vje, mje, cjc, vjc, mjc, fc, tf, tr, xtf, itf, vtf;
        double rb, rbm, irb;
    };
    void bjtDevice(int c, int b, int e, int externalBase, int sign, const BJT& m) {
        const Jet be{sign * (voltage(b) - voltage(e)), 1, 0};
        const Jet bc{sign * (voltage(b) - voltage(c)), 0, 1};
        const auto f = junction(be, m.nf * thermal, m.is);
        const auto r = junction(bc, m.nr * thermal, m.is);
        auto invBase = Jet{1} - be * (m.var ? 1 / m.var : 0)
                              - bc * (m.vaf ? 1 / m.vaf : 0);
        if (m.ikf || m.ikr) {
            const auto highInjection = Jet{1} + (f * (m.ikf ? 1 / m.ikf : 0)
                                              + r * (m.ikr ? 1 / m.ikr : 0)) * 4;
            const auto correction = highInjection.value <= 0 ? Jet{1}
                : m.nk == .5 ? squareRoot(highInjection) : power(highInjection, m.nk);
            invBase = invBase / (Jet{.5} + correction * .5);
        }
        const auto t = (f - r) * invBase;
        auto ibf = f * (1 / m.bf) + junction(be, m.ne * thermal, m.ise) + be * 1e-12;
        auto ibr = r * (1 / m.br) + junction(bc, m.nc * thermal, m.isc) + bc * 1e-12;
        if (b != externalBase) {
            double rb = m.rbm + (m.rb - m.rbm) * invBase.value;
            if (m.irb > 0) {
                const double current = std::max((ibf.value + ibr.value) / m.irb, 1e-9);
                const double root = std::sqrt(1 + 14.59025 * current);
                const double angle = 14.59025 * std::sqrt(current) / (2.4317 * (root + 1));
                const double tangent = std::tan(angle);
                rb = m.rbm + 3 * (m.rb - m.rbm) * (tangent - angle) / (angle * tangent * tangent);
            }
            conductance(externalBase, b, 1 / rb, true);
        }
        if (m.cje || m.cjc || m.tf || m.tr) {
            auto qf = depletion(be, m.cje, m.vje, m.mje, m.fc);
            const auto qr = depletion(bc, m.cjc, m.vjc, m.mjc, m.fc) + r * m.tr;
            if (m.tf != 0) {
                auto forward = f;
                if (be.value > 0) {
                    const double k = m.vtf > 0 ? 1 / (1.44 * m.vtf) : 0;
                    const double expV = m.xtf * std::exp(std::clamp(bc.value * k, -70.0, 70.0));
                    auto bias = Jet{expV, 0, expV * k};
                    if (m.itf > 0) {
                        const auto current = f / (f + Jet{m.itf});
                        bias = bias * current * current;
                    }
                    forward = f * (Jet{1} + bias) * invBase;
                }
                qf = qf + forward * m.tf;
            }
            ibf = ibf + Jet{companion(qf.value), inverse * qf.p, inverse * qf.q};
            ibr = ibr + Jet{companion(qr.value), inverse * qr.p, inverse * qr.q};
        }
        terminal(c, c, b, e, sign, t - ibr);
        terminal(b, c, b, e, sign, ibf + ibr);
        terminal(e, c, b, e, sign, (t + ibf) * -1);
    }

    std::vector<Companion> companions;
    std::size_t nextCompanion = 0;
    double inverse = 0;
    int unknowns;
};
}
