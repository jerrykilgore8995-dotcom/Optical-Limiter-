* Reduced sidechain. The full vendor deck is the measurement reference.
Bselect sc 0 V={FF*V(pre)+(1-FF)*V(gr)}
Escbuf scbuf 0 sc 0 1
CH scbuf hp 100n
RH hp 0 {1/(6.283185307*max(HPF,10)*100n)}
Bhpsel hpsel 0 V={HPF>0 ? V(hp) : V(scbuf)}
Ehfin hfin 0 hpsel 0 -1
RHF hfin hfout 10k
CFB hfout 0 220p
* Precision full-wave rectifier. Both half cycles drive the original detector loads.
Bhalfin halfin 0 V={max(V(hfout),0)}
RHALF halfin half 12k
CHB half 0 22p
Brectin rectin 0 V={2*V(half)-V(hfout)}
RRECT rectin rectpole 9.98k
CRECT rectpole 0 22p
ERECT rect 0 rectpole 0 1
* Average and peak detector branches retain their independent history.
RA rect avg 10k
CA avg 0 100n
Eavg avgin 0 avg 0 1.571
RAF avgin avgo 5.71k
CAF avgo 0 22p
* The closed precision loop reduces incremental charging resistance.
* Keep the 100-ohm current limit and the measured peak-release network.
* Compact two-state fit to peak-amplifier recovery at the rectifier transition.
* The band-pass and high-pass terms reproduce the positive charging pulse.
CPKFB pkcomp 0 680p
Cpkvel pkvel 0 680p
Gpkint 0 pkcomp pkvel 0 0.0001875708347
Gpkfeed 0 pkvel hfout 0 0.0001875708347
Gpkrestore pkvel 0 pkcomp 0 0.0001875708347
Gpkdamp pkvel 0 pkvel 0 0.0006098911232
Epkedge pkedge 0 POLY(3) (hfout,0) (pkcomp,0) (pkvel,0) 0 8.69355234 -8.69355234 -22.793954
* Fit the edge strength over the measured range, including weak pulses and compression at high drive.
Bcharge 0 hold I={min(max(max(V(rect),1.02*V(pkedge)+max(0,min(0.45*(V(pkedge)-0.04),0.28*(0.18-V(pkedge))))-0.57*max(V(pkedge)-0.62,0))-V(hold),0)/8,max(13.5-V(hold),0)/100)}
CP hold 0 100n
RP hold 0 47k
Epkbuf pkout 0 hold 0 1
Bdet det 0 V={RATIO<8 ? max(V(avgo),1e-9) : max(V(pkout),1e-9)}
VPOT pw 0 {3*max(0.01,50000*SENS/30)/(max(0.01,50000*(1-SENS/30))+max(0.01,50000*SENS/30))}
* Below 1 mV the detector is far below either knee at every threshold.
* A bounded log derivative there avoids false Newton retries at silence.
Bsub subin 0 V={min(14,max(-14,-2*log10(max(1m,min(15,V(det))))-V(pw)))}
RSUB subin subout 10k
CSUB subout 0 22p
Einv invin 0 subout 0 -1
RIF invin ov100 10k
CIF ov100 0 22p
Bover over 0 V={10*max(-15,min(15,V(ov100)))}
Bknee knee 0 V={SENS<=0 ? 0 : (FF>0.5 ? SOFT : 0) ? (V(over)<-5 ? 0 : V(over)>5 ? V(over) : pow(V(over)+5,2)/20) : max(V(over),0)}
BRR rrin 0 V={-0.1*V(knee)*(FF ? 1-1/RATIO : RATIO-1)}
RRF rrin rrout {10000*(FF ? 1-1/RATIO : RATIO-1)}
CRRF rrout 0 22p
Bcmd cmd 0 V={min(4,max(0,-V(rrout)))}
CT timing 0 1u
RTIME cmd timing {REL/1u}
BT cmd timing I={max(V(cmd,timing),0)*(1u/ATT-1u/REL)}
Etc tc 0 timing 0 1
Breq req 0 V={cg((exp(1.151292546497023*min(4,max(0,V(tc))))-1)/10000)}
Bdrive drive 0 V={min(3,max(0,ternary_fcn(cg(V(req))<5e-5,0.1*sqrt(max(1e-20,(cg(V(req))-2e-8)/(5e-5-2e-8))),ternary_fcn(cg(V(req))<1/600,0.1*pow(cg(V(req))/5e-5,1/1.52288),pow(cg(V(req))*600,1/0.79248)))))}
Esense sense 0 drive 0 1
