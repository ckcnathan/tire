
clear AMBTMP datamode ET fileid FX FY FZ IA MX MZ N NFX NFY P pathname RE RL RST SA SL SR testid tireid TSTC TSTI TSTO V
clc


%datamode = LongOrLat(); % Choose longetudinal or lateral mode
[fileid, pathname] = uigetfile({'*.mat;*.dat'},'File Selector');
if( isa(fileid, 'double') && isequal(fileid, 0) )
    % No file was chosen, so we return
    return
end

load([pathname fileid]); %load METRIC .mat file only (for now)

%%
close
nexttile
plot(ET,P) 
nexttile    
plot(ET,IA)
nexttile
plot(ET,FZ)
nexttile
hold on
plot(ET,SA)

%% remove stuff
x = 19931;
x1 = 52322;

ET(x:x1) = [];
FY(x:x1) = [];
FX(x:x1) = [];
FZ(x:x1) = [];
IA(x:x1) = [];
MX(x:x1) = [];
MZ(x:x1) = [];
N(x:x1) = [];
NFX(x:x1) = [];
NFY(x:x1) = [];
P(x:x1) = [];
RE(x:x1) = [];
RL(x:x1) = [];
RST(x:x1) = [];
SA(x:x1) = [];
SL(x:x1) = [];
SR(x:x1) = [];
TSTC(x:x1) = [];
TSTI(x:x1) = [];
TSTO(x:x1) = [];
V(x:x1) = [];

%% cat stuff
ET1 = cat(1, ET1, ET);
FY1 = cat(1, FY1, FY);
FX1 = cat(1, FX1, FX);
FZ1 = cat(1, FZ1, FZ);
IA1 = cat(1, IA1, IA);
MX1 = cat(1, MX1, MX);
MZ1 = cat(1, MZ1, MZ);
N1 = cat(1, N1, N);
NFX1 = cat(1, NFX1, NFX);
NFY1 = cat(1, NFY1, NFY);
P1 = cat(1, P1, P);
RE1 = cat(1, RE1, RE);
RL1 = cat(1, RL1, RL);
RST1 = cat(1, RST1, RST);
SA1 = cat(1, SA1, SA);
SL1 = cat(1, SL1, SL);
SR1 = cat(1, SR1, SR);
TSTC1 = cat(1, TSTC1, TSTC);
TSTI1 = cat(1, TSTI1, TSTI);
TSTO1 = cat(1, TSTO1, TSTO);
V1 = cat(1,V1,V);


%% save stuff first 
ET1 = ET;
FY1 = FY;
FX1 = FX;
FZ1 = FZ;
IA1 = IA;
MX1 = MX;
MZ1 = MZ;
N1 = N;
NFX1 = NFX;
NFY1 = NFY;
P1 = P;
RE1 = RE;
RL1 = RL;
RST1 = RST;
SA1 = SA;
SL1 = SL;
SR1 = SR;
TSTC1 = TSTC;
TSTI1 = TSTI;
TSTO1 = TSTO;
V1 = V;

%% before trimming
clear ET FY FX FZ IA MX MZ N NFX NFY P RE RL RST SA SL SR TSTC TSTI TSTO V
ET = ET1;
FY = FY1;
FX = FX1;
FZ = FZ1;
IA = IA1;
MX = MX1;
MZ = MZ1;
N = N1;
NFX = NFX1;
NFY = NFY1;
P = P1;
RE = RE1;
RL = RL1;
RST = RST1;
SA = SA1;
SL = SL1;
SR = SR1;
TSTC = TSTC1;
TSTI = TSTI1;
TSTO = TSTO1;
V = V1;

%% save data
clear ET1 FY1 FX1 FZ1 IA1 MX1 MZ1 N1 NFX1 NFY1 P1 RE1 RL1 RST1 SA1 SL1 SR1 TSTC1 TSTI1 TSTO1 V1
clear x x1 datamode 
save HOOSIER_43075_20.5x7.5-13_R20_7rim_trimmed_long.mat
% save GOODYEAR_D2704_20.0x7.2-13_7rim_long_trimmed.mat

%% (117292:124223)

close
nexttile
plot(FZ)
nexttile
plot(P)
nexttile
plot(FY)
nexttile
plot(IA)
nexttile
plot(SA)
