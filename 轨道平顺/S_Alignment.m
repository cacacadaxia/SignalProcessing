function Sv=S_Alignment(Omega)
k=0.25;
Aa=0.0339*100;
Omega_c=0.8245;
Sv=(k.*Aa.*Omega_c.^2)./(Omega.^2.*(Omega.^2+Omega_c.^2));
