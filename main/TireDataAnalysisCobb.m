%%
% * AT: Ambient Temperature [degC]
% * B,C,D,E: Magic Formula Coefficients
% * ET: Elapsed Time [s]
% * FX: Traction Force [N]
% * FY: Lateral Force [N]
% * FZ: Vertical Force [N]
% * IA: Inclination Angle [deg]
% * MX: Overturning Moment [Nm]
% * MZ: Aligning Moment [Nm]
% * N: Wheel rotational speed [rpm]
% * NFX: Normal Force in X [normalized against vectorial load]
% * NFY: Normal Force in Y [normalized against vectorial load]
% * P: Pressure [kPa]
% * RE: Effective Radius [cm]
% * RL: Loaded Radius [cm]
% * RST: Road Surface Temperature [degC]
% * SA: Slip Angle [deg]
% * SR: Slip Ratio [-]
% * TSTC: Tire Surface Temperature Center [degC]
% * TSTI: Tire Surface Temperature Inner [degC]
% * TSTO: Tire Surface Temperature Outer [degC]
% * V: Forward Speed of Axle [m/s]
%% Initialization
% clear
clearvars -except tire;
clc;
close all;

addpath(fullfile(pwd,'Functions'));

datamode = LongOrLat(); % Choose longetudinal or lateral mode
[fileid, pathname] = uigetfile({'*.mat;*.dat'},'File Selector');
if( isa(fileid, 'double') && isequal(fileid, 0) )
    % No file was chosen, so we return
    return
end

[AT, ET, FX, FY, FZ, IA, MX, MZ, N, NFX, NFY, P, RE, RL, RST, SA, SR, TSTC, TSTI, TSTO, V]...
    = ImportRawData( fullfile(pathname, fileid) );
%%
% Functions for converting between lbf and N

lbf2N = @(lbf)lbf*4.4482216152605;
N2lbf = @(N)N/4.4482216152605;
%%
% Functions for converting between kPa and psi

psi2kPa = @(psi)psi/0.145037737730217;
kPa2psi = @(kPa)kPa*0.145037737730217;
%%
%
%% Identify data points tested
% Tires are tested following pre-defined test routines. To calculate the Magic
% Formula coefficients it is important to identify test points of *normal force,
% tire pressure and inclination angle*. The correct test values can be identified
% in the following figures:
%
% *Normal force*:
%
% You can see in the example below that the tire was tested under different
% normal forces. The Magic Formula parameters depend on the normal force FZ of
% the tire, where the normal force should be in Newtons. The values are around
% 50,100,150,250 pounds (1pound = 4.448N). We save those values to create our
% normal force bins later.
%
% *Figure to identify the Variation of FZ:*

figure
plot(ET,FZ,'.');
grid on
xlabel('Elapsed Time (s)');
ylabel('FZ (N)');
title('Variation of FZ','FontSize',10);
[countsFZ,edgesFZ] = histcounts(N2lbf(FZ));
countsFZ = movmean(countsFZ,5);
countsFZ(countsFZ<100) = 0;
[~,locsFZ] = findpeaks([countsFZ(2), countsFZ, countsFZ(end-1)]);
%%
% *Tire pressure:*
%
% Tire pressure is often denoted in kPa but to visualize the data it was converted
% to pound-force per square inch (1 kPa = 0.145 psi).

figure
plot(ET,0.145*P);
xlabel('Elapsed Time (s)');
ylabel('P (psi)');
grid on
title('Variation of Pressure','FontSize',10);
[countsP,edgesP] = histcounts(kPa2psi(P));
countsP = movmean(countsP,3);
[~,locsP] = findpeaks([countsP(2), countsP, countsP(end-1)]);
%%
% *Inclination Angle:*
%
% The inclination angle is the angle of the wheel relative to the surface.

figure
plot(ET,IA);
xlabel('Elapsed Time (s)');
ylabel('IA ({\circ}''s)');
grid on
title('Variation of IA','FontSize',10);
[countsIA,edgesIA] = histcounts(IA);
% countsIA = movmean(countsIA,3);
countsIA(countsIA<100) = 0;
[~,locsIA] = findpeaks([countsIA(2), countsIA, countsIA(end-1)]);

%% Slip angle for combined

[countsSA, edgesSA] = histcounts(SA);
[~,locsSA] = findpeaks([countsSA(2), countsSA, countsSA(end-1)]);

%%
% *Note:* If you want to check the data that was imported is proper comment
% out the code after this point and check the graphs above
%% Create bins for data
% To group data into bins, logic statements are used. Those bins describe which
% values of normal force, pressure or inclination angle have been tested. Tolerances
% are applied to accomodate measured data into bins.
%
% *Note:* The choice of test values has to be done manually.

FZ_binvalues = lbf2N(unique(round(abs(edgesFZ(locsFZ)))));
P_binvalues = unique(round(edgesP(locsP)));
IA_binvalues = unique(round(edgesIA(locsIA)));
SA_binvalues = unique(round(edgesSA(locsSA)));

if( isempty(FZ_binvalues) || isempty(P_binvalues) || isempty(IA_binvalues) )
    error('One or more of the *_binValues arrays was empty.')
end
%%
% *Pressure bins*

P_eps = 0.8;  % Pressure tolerance (psi)
P_bin = kPa2psi(P)>(P_binvalues-P_eps-0.3)&kPa2psi(P)<(P_binvalues+P_eps); %(psi)
%%
% *Note*: Using implicit expansion (available from release 2016b on). Suggested
% to use repmat() when working with older releases.
%
% *Inclination angle bin*

IA_eps = 0.2; % Inclination angle tolerance (deg)
IA_bin = (IA>IA_binvalues-IA_eps)&(IA<IA_binvalues+IA_eps); %(deg)
%%
% *Normal force bin*

FZ_eps1 = lbf2N(35); % Normal force tolerance (N)
FZ_bin = abs(FZ)>((FZ_binvalues-FZ_eps1))&abs(FZ)<((FZ_binvalues+FZ_eps1)); %(N)

%% Slip angle bin

SA_eps = 0.1;
SA_bin = SA>((SA_binvalues-SA_eps))&SA<((SA_binvalues+SA_eps));
%%
% *Zero slip bin*
%
% The zero slip bin is used where the slip values should be close to zero but
% these values depends on whether the model is longitudinal or lateral. If the
% data is longitudinal slip angles are used and slip ratio for lateral.

if datamode == 1
    S_0 = (-1<SA)&(SA<1);
else
    S_0 = (-1<SR)&(SR<1);
end

%% Data analysis
% To obtain Magic Formula coefficients for varying data points the code loops
% through tire pressure, normal force and inclination angle bins.
%
% *Loop through pressure bins*

q=0;
for i=1:length(P_binvalues)
    for m=1:size(FZ_bin,2)
        for n=1:size(IA_bin,2)
            validIdx = FZ_bin(:,m) & P_bin(:,i) & IA_bin(:,n);
            ET_binfzia{i,n,m}  =  ET(validIdx);  % Time Bins
            FZ_binfzia{i,n,m}  =  FZ(validIdx);  % Vertical Load bins
            IA_binfzia{i,n,m}  =  IA(validIdx);  % Inclination Angle bins
            MX_binfzia{i,n,m}  =  MX(validIdx);  % Overturning Moment bins
            MZ_binfzia{i,n,m}  =  MZ(validIdx);  % Aligning Moment bins

            if datamode == 1
                F_binfzia{i,n,m}  =  -FX(validIdx);  % Logitudinal Force Bins
                NF_binfzia{i,n,m} =  NFX(validIdx); % Force Coefficient bins
                S_binfzia{i,n,m}  =  -SR(validIdx);  % Slip Ratio Bins
            elseif datamode == 2
                F_binfzia{i,n,m}  =  FY(validIdx);  % Lateral Force Bins
                NF_binfzia{i,n,m} =  NFY(validIdx); % Force Coefficient bins
                S_binfzia{i,n,m}  =  -SA(validIdx);  % Slip Angle Bins
            end

            if datamode ==2
                sp_f{i,n,m}=csaps(S_binfzia{i,n,m},F_binfzia{i,n,m},.1);
                sp_mz{i,n,m}=csaps(S_binfzia{i,n,m},MZ_binfzia{i,n,m},.1);
                sp_mx{i,n,m}=csaps(S_binfzia{i,n,m},MX_binfzia{i,n,m},.1);
                %sp_rl{i,n,m}=csaps(S_binfzia{i,n,m},rl,.1);
            elseif datamode ==1
                sp_f{i,n,m}=csaps(S_binfzia{i,n,m},F_binfzia{i,n,m},.995);
                S_SH = mean(fnzeros(sp_f{i,n,m}));
                S_binfzia{i,n,m} = S_binfzia{i,n,m}-S_SH;
                sp_f{i,n,m}=csaps(S_binfzia{i,n,m},F_binfzia{i,n,m},.995);

            end

            % figure('Name',[ ': Aligning Moment vs. Slip Angle & Vertical Load'], 'numbertitle','off')
            % %subplot(3,1,1)
            % hold on
            % plot(S_binfzia{i,n,m},F_binfzia{i,n,m},'.','color',[1 .5 0])
            % fnplt(sp_f{i,n,m},'b')
            % title({['Fz= ' num2str(round(mean(FZ_binvalues(m)))) ' N'];['IA= ' num2str(round(mean(IA_binvalues(n)))) '°']})
            % xlabel('Slip Angle')
            % ylabel('Lateral Force')
            % line([min(S_binfzia{i,n,m}) max(S_binfzia{i,n,m})],[0 0],'color','k')
            % line([0 0],[min(F_binfzia{i,n,m}) max(F_binfzia{i,n,m})],'color','k')
            % legend('Test Data','Fitted Data','location','southeast')
            % close
            %subplot(3,1,2)
            % hold on
            % plot(S_binfzia{i,n,m},MZ_binfzia{i,n,m},'.','color',[.5 .5 .5])
            % fnplt(sp_mz{i,n,m},'b')
            % xlabel('Slip Angle')
            % ylabel('Aligning Moment')
            % line([min(S_binfzia{i,n,m}) max(S_binfzia{i,n,m})],[0 0],'color','k')
            % line([0 0],[min(MZ_binfzia{i,n,m}) max(MZ_binfzia{i,n,m})],'color','k')
            % subplot(3,1,3)
            % hold on
            % plot(S_binfzia{i,n,m},MX_binfzia{i,n,m},'.','color',[.5 .5 .5])
            % fnplt(sp_mx{i,n,m},'b')
            % xlabel('Slip Angle')
            % ylabel('Overturning Moment')
            % line([min(S_binfzia{i,n,m}) max(S_binfzia{i,n,m})],[0 0],'color','k')
            % line([0 0],[min(MX_binfzia{i,n,m}) max(MX_binfzia{i,n,m})],'color','k')

            if datamode == 2
                for xs = -13:1:13
                    q=q+1;
                    fmdata(q,1)=xs;
                    fmdata(q,2)=IA_binvalues(n);
                    fmdata(q,3)=FZ_binvalues(m);
                    fmdata(q,4)=fnval(sp_f{i,n,m},xs);
                    fmdata(q,5)=fnval(sp_mz{i,n,m},xs);
                    fmdata(q,6)=fnval(sp_mx{i,n,m},xs);
                end
            elseif datamode ==1
                for xs = -0.21:0.014:0.21
                    q=q+1;
                    fmdata(q,1)=xs;
                    fmdata(q,2)=IA_binvalues(n);
                    fmdata(q,3)=FZ_binvalues(m);
                    fmdata(q,4)=fnval(sp_f{i,n,m},xs);
                end
            end

            % if n==1
            % % plot3(S_binfzia{i,n,m},FZ_binvalues(m)*ones(length(S_binfzia{i,n,m})),F_binfzia{i,n,m})
            % plot3(fmdata(q,1),FZ_binvalues(m)*ones(length(fmdata(q,1))),fmdata(q,4))
            % hold on

        end
    end


    fmdata = sortrows(fmdata,[2,1,3]);
    incls = unique(round(fmdata(:,2)))' ;
    nincls = length(incls);
    % slips = unique(round(fmdata(:,1)))' ;
    slips = unique(fmdata(:,1))' ;
    nslips = length(slips);

    for nn = 1:length(IA_binvalues)

        inx0 = find(fmdata(:,2) == IA_binvalues(nn)); % zero camber points
        % then grab the rest of the zero camber condition data:
        fmdata0 = fmdata(inx0,:);
        % Next, we transpose the arrays to make our spline functions happy:
        loads = mean(reshape(fmdata0(:,3),[],nslips),2)';
        nloads = length(loads);
        %Take a look at FZ:
        fz0 = reshape(fmdata0(:,3),nloads,nslips)';
        fy0 = reshape(fmdata0(:,4),nloads,nslips)';
        if datamode==2
            mz0 = reshape(fmdata0(:,5),nloads,nslips)';
            mx0 = reshape(fmdata0(:,6),nloads,nslips)';
        end
        nfy0 = fy0./fz0;

        %%
        if datamode ==2
            tire.lat.FY{i,nn} = csaps({slips,loads},fy0);
            tire.lat.MZ{i,nn} = csaps({slips,loads},mz0);
            tire.lat.MX{i,nn} = csaps({slips,loads},mx0);
            tire.lat.FZval = FZ_binvalues;
            tire.lat.Pval = P_binvalues;
            tire.lat.IAval = IA_binvalues;
        elseif datamode ==1
            tire.long.FX{i,nn} = csaps({slips,loads},fy0);
            tire.long.FZval = FZ_binvalues;
            tire.long.Pval = P_binvalues;
            tire.long.IAval = IA_binvalues;
            % figure('Name',[': Lateral Force vs. Slip Angle & Vertical Load '])
            % fnplt(LATE_SLIP_VERT)
            % alpha 0.2
            % xlabel('Slip Angle (deg)')
            % ylabel('Vertical Load (N)')
            % zlabel('Lateral Force (N)')
            % view(45,45)
            % hold on
        end
    end

    clear fmdata
    q=0;
end

tire.name = 'HOOSIER_20.5x7.5-13_R20_7rim';
%tire.name = 'GOODYEAR 20.0x7.0-13 R065 7"rim';
tire.lat.legend = 'pressure x IA';
tire.long.legend = 'pressure x IA';


%%
spline = tire;
% save HOOSIER_20.5x7.5-13_R20_7rim.mat spline
% save GOODYEAR_20.0x7.0-13_R065_7rim.mat spline

%%
% clc
% i = 2;
% n = 1;
% m = 5;
% plot(S_binfzia{i,n,m},F_binfzia{i,n,m},'.','color',[.5 .5 .5])
% hold on
% xss = 0:0.1:13;
% for ff = 1:length(xss)
%     xff(ff) = fnval(tire.lat.FY{i,n},[xss(ff);1500]);
% end
% plot(xss,xff)

%%
clc
clear xfy xs sh
i = 1;  %pressure   
n = 1;  %IA
m = 4;  %F

plot(S_binfzia{i,n,m},F_binfzia{i,n,m},'.','Color',[0, 0, 0.8, 0.2])
hold on

xs = 0:0.1:12;
sh = 0;
for k=1:length(xs)
    xfy(k) = fnval(tire.lat.FY{i,n},[xs(k);FZ_binvalues(m)]);
    % if xs(k) == 0
    %     sh = xfy(k);
    %     xfy(k) = xfy(k)-sh;
    % end
end
plot(xs,xfy,'LineWidth',2)
axis([0 13 0 2000])
grid on
legend("raw test data","fitted spline")
title("Spline Tire Model Fit")
ylabel("FY [N]")
xlabel("Slip Angle [deg]")

% clear xss xff
% xss = 0:0.01:0.2;
% for ff = 1:length(xss)
%     xff(ff) = fnval(tire.long.FX{i,n},[xss(ff);FZ_binvalues(m)]);
% end
% plot(xss,xff)
