# tire
For testing and fitting tire data from the Tire Testing Consortium (TTC). Script written for Concordia Formula Racing 2023-2024. <br/> 

Data parser is based on https://blogs.mathworks.com/student-lounge/2018/08/01/analyzing-tire-data/  <br/>
Tire data is fitted to a cubic spline, partially based on Bill Cobb's method from the TTC forum.  <br/>

TireDataAnalysisCobb.m -> main tire fitting script.  <br/>
manualtrimmer.m           -> used for trimming TTC data so the parser in TireDatAnalysisBillyBob.m can work properly.  <br/>
tireplots.m               -> for plotting and comparing results.  <br/>              
