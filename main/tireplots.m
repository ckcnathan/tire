
clearvars -except spline
clc

FZ = [1100];   % 709, 323.87, 768, 425.98

for i=1:length(FZ)
    q=1;
    sh = 1*fnval(spline.lat.FY{1,1},[0;FZ(i)]);
    for ss = 0:0.5:13
        xss(i,q) = ss;
        ffy(i,q) = (1*fnval(spline.lat.FY{1,1},[ss;FZ(i)])-sh)/FZ(i);
        q=q+1;
    end
    plot(xss(i,:),ffy(i,:),"color","#00ff57")
    hold on
end
%%
grid on

title("Lateral mu vs Slip Angle at Various Normal Loads, 0 camber, 10psi")
% legend("10inch, Fz = 300N","10inch, Fz = 700N","10inch, Fz = 1100N","13inch, Fz = 300N","13inch, Fz = 700N","13inch, Fz = 1100N","location","southeast")
legend("10inch","13inch","location","southeast")

xlabel("Slip angle (deg)")
ylabel("Lateral mu")
%%
% plot(xss,(ffy(1,:)+ffy(3,:))/2)
% hold on
% plot(xss,(ffy(2,:)*2)/2)
% grid on
% 
% title("Effective Axle lateral mu with 1 lateral g, 56% load transfer distribution")
% xlabel("Slip Angle [deg]")
% ylabel("Effective Axle lateral mu")
% legend("With Load Transfer","Without Load Transfer","location","southeast")