clc
clear
format long g;

Lambda=1:0.5:30;
Omega=2.*pi./Lambda;
Sv=S_Alignment(Omega);
loglog(Lambda,Sv),shg
xlim([0 30]);
ylim([1e-4 1e3]);
set(gca,'XDir','reverse');
grid on;

function Sv=S_Alignment(Omega)
k=0.25;
Aa=0.0339*100;
Omega_c=0.8245;
Sv=(k.*Aa.*Omega_c.^2)./(Omega.^2.*(Omega.^2+Omega_c.^2));
end