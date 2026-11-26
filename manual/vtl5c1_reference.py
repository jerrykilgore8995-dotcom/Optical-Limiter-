"""Datasheet-anchor fit only; not independent validation or measured THD."""
import math
DARK=50e6
def conductance(ma):
    ma=max(0.,ma)
    if ma<1:
        return 1/DARK+(1/20000-1/DARK)*ma**2
    if ma<10:
        return (1/20000)*ma**(math.log(20000/600)/math.log(10))
    return min(1/200,(1/600)*(ma/10)**(math.log(600/200)/math.log(4)))
def gr(g):
    return 20*math.log10(1+10000*g)
# Exponential conductance surrogate calibrated to 100k endpoint at 35ms
# from 40mA. Datasheet turn-off is a maximum, not a typical time constant.
tau=0.035/math.log((1/200-1/DARK)/(1/100000-1/DARK))
def release(g,t):
    return 1/DARK+(g-1/DARK)*math.exp(-t/tau)
if __name__=="__main__":
    for ma,expected in [(0,DARK),(1,20000),(10,600),(40,200)]:
        g=conductance(ma); r=1/g
        assert abs(r/expected-1)<1e-9
        print(f"{ma:2} mA R={r:.3f} ohm GR={gr(g):.4f} dB")
    print(f"Endpoint-fit tau={tau*1000:.3f} ms")
    print(f"40mA release R at 35ms={1/release(conductance(40),.035):.1f} ohm")
    # Do not claim distortion: voltage-dependent conductance data unavailable.
