#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cmath>
namespace optical {
constexpr double thermal=1.380649e-23*300.15/1.602176634e-19;
struct Jet {
    double value, p = 0.0, q = 0.0;
};
inline Jet operator+(Jet a, Jet b) { return {a.value+b.value,a.p+b.p,a.q+b.q}; }
inline Jet operator-(Jet a, Jet b) { return {a.value-b.value,a.p-b.p,a.q-b.q}; }
inline Jet operator*(Jet a, double b) { return {a.value*b,a.p*b,a.q*b}; }
inline Jet operator*(Jet a, Jet b) {
    return {a.value*b.value,a.p*b.value+a.value*b.p,a.q*b.value+a.value*b.q};
}
inline Jet operator/(Jet a, Jet b) {
    const double inv = 1.0/b.value, y = a.value*inv;
    return {y,(a.p-y*b.p)*inv,(a.q-y*b.q)*inv};
}
inline Jet squareRoot(Jet x) {
    const double value = std::sqrt(x.value), slope = 0.5/value;
    return {value,x.p*slope,x.q*slope};
}

inline double junctionExp(double x) {
    // The forward-junction load bounds x to [-3, 70]. Range reduction
    // keeps r within half a log(2), where this degree-eight series is
    // accurate to better than one part in a billion.
    const int exponent=static_cast<int>(x*1.4426950408889634+4.5)-4;
    const double r=x-exponent*0.6931471805599453;
    const double square=r*r, fourth=square*square;
    const double low=1+r+square*(0.5+r/6);
    const double high=(1.0/24+r/120)+square*(1.0/720+r/5040);
    const auto power=std::bit_cast<double>(static_cast<std::uint64_t>(exponent+1023)<<52);
    return (low+fourth*(high+fourth/40320))*power;
}

inline Jet junction(Jet volts, double emissionVoltage, double saturation) {
    if (saturation == 0) return {0};
    const double x = volts.value/emissionVoltage;
    if (x < -3) {
        const double a = 3.0/(2.718281828459045*x);
        const double cube = a*a*a, slope = saturation*3*cube/volts.value;
        return {-saturation*(1 + cube), slope*volts.p, slope*volts.q};
    }
    // The linear continuation above 1.8 V keeps a rejected Newton
    // iterate finite. It is outside the normal operating region.
    const double limited = std::min(x, 70.0);
    const double exp = junctionExp(limited), slope = saturation*exp/emissionVoltage;
    return {saturation*(exp*(1 + x-limited)-1), slope*volts.p, slope*volts.q};
}

template<double exponent>
const auto junctionPowerTable=[] {
    std::array<double,8193> data{};
    for (std::size_t i=0;i<data.size();++i) data[i]=std::pow(0.5+double(i)/256,-exponent);
    return data;
}();

template<double exponent>
inline Jet capacitance(Jet volts, double zero, double potential, double fc) {
    if (zero == 0) return {0};
    const double knee = potential*fc;
    double value, slope;
    if (volts.value < knee) {
        const double x = 1-volts.value/potential;
        const double index=(x-0.5)*256;
        double m;
        if (index>=0 && index<8192) {
            const int slot=static_cast<int>(index);
            const auto& data=junctionPowerTable<exponent>;
            m=data[slot]+(index-slot)*(data[slot+1]-data[slot]);
        } else m=std::pow(x,-exponent);
        value = potential*zero*(1-x*m)/(1-exponent);
        slope = zero*m;
    } else {
        const double f = std::pow(1-fc,-exponent);
        const double offset = volts.value-knee;
        const double curvature = f*exponent/(potential*(1-fc));
        value = potential*zero*(1-(1-fc)*f)/(1-exponent)
                + zero*offset*(f+0.5*offset*curvature);
        slope = zero*(f+offset*curvature);
    }
    return {value,slope*volts.p,slope*volts.q};
}


}
