# tire
For testing and fitting tire data from the Tire Testing Consortium (TTC).

Data parser is based on https://blogs.mathworks.com/student-lounge/2018/08/01/analyzing-tire-data/
tire data is fitted to a to a cubic spline, partially based on Bill Cobb's method from the TTC forum.

TireDatAnalysisBillyBob.m -> main tire fitting script
manualtrimmer.m           -> used for trimming TTC data so the parser in TireDatAnalysisBillyBob.m can work properly.
tireplots.m               -> for plotting and comparing results            
