    % 
    % val(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2
    % x = IA_binvalues;
    % y = FZ_binvalues;
    % P_binvalues
    %
    % S_H_surf_IA_P{i} 
    % S_V_surf_IA_P{i}
    % Mu_surf_IA_P{i} 
    % CS_surf_IA_P{i}  
    % B_surf_IA_P{i}   
    % C_surf_IA_P{i} 
    % D_surf_IA_P{i} 
    % E_surf_IA_P{i}
    clc
    close
    datamode = 2;
    FZ = 700;
    IA = 0;
    
    CS_bar = MF.lat.CS_surf{1};
    S_H_bar = MF.lat.S_H_surf{1};
    mu_bar = MF.lat.Mu_surf{1};

    if datamode == 2
        Slip = -12:0.1:12;
        Slip_bar = (CS_bar)./mu_bar./FZ.*(tan((Slip-S_H_bar)*pi/180));
    elseif datamode == 1
        Slip = -0.25:0.01:0.25;
        Slip_bar = (CS_bar.*Slip-S_H_bar)./mu_bar./FZ;
    end

    F_bar = MagicFormula([B_surf_IA_P{2}(IA,FZ),E_surf_IA_P{2}(IA,FZ)],Slip_bar);
    F = F_bar.*mu_bar.*FZ;
    % F = 1/mu_bar.*(F_bar./FZ);

    plot(Slip,F)
    




