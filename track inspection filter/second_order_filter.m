%二阶补偿滤波器 级联 抗混滤波器

%抗混滤波器
k = 10^5 / 2^14;
a = k/2;
b = (k^2-a^2)^(1/2);

v   = 16;
pesi= 0.0001 : 0.0001 : 10;
Fs1 = (a^2+b^2) ./(-(2*pi*pesi*v).^2 + 2*1j*a*2*pi.*pesi*v + (a^2+b^2));
v   = 72;
Fs2 = (a^2+b^2) ./(-(2*pi*pesi*v).^2 + 2*1j*a*2*pi.*pesi*v + (a^2+b^2));
v   = 128;
Fs3 = (a^2+b^2) ./(-(2*pi*pesi*v).^2 + 2*1j*a*2*pi.*pesi*v + (a^2+b^2));

%% 补偿滤波器
%% 滤波器系数
w = 0.0001 : 0.0001 :10;
k = 10^5 / 2^14;
a = k/2;
b = (k^2-a^2)^(1/2);
v = 16;
t = 0.25/v;
a1 = ( 1 + a*t + a^2*t*t )/( t^2*(a^2+b^2) );
a2 = 2*(1-a^2*t^2) / ( t^2*(a^2+b^2) ) ;
a3 = ( 1 - a*t + a^2*t*t ) / ( t^2*(a^2+b^2) ) ;
Djw1 = a1-a2*exp(-1j*2*pi*0.25.*w)+a3*exp(-2j*2*pi*0.25.*w);

v  = 72;
t  = 0.25/v;
a1 = ( 1 + a*t + a^2*t*t )/( t^2*(a^2+b^2) );
a2 = 2*(1-a^2*t^2) / ( t^2*(a^2+b^2) ) ;
a3 = ( 1 - a*t + a^2*t*t ) / ( t^2*(a^2+b^2) ) ;
Djw2 = a1-a2*exp(-1j*2*pi*0.25.*w)+a3*exp(-2j*2*pi*0.25.*w);

v  = 128;
t  = 0.25/v;
a1 = ( 1 + a*t + a^2*t*t )/( t^2*(a^2+b^2) );
a2 = 2*(1-a^2*t^2) / ( t^2*(a^2+b^2) ) ;
a3 = ( 1 - a*t + a^2*t*t ) / ( t^2*(a^2+b^2) ) ;
Djw3 = a1-a2*exp(-1j*2*pi*0.25.*w)+a3*exp(-2j*2*pi*0.25.*w);

%% 
figure
semilogx(w,20*log10(Fs1),w,20*log10(Fs2),w,20*log10(Fs3));
ylabel('dB ');
xlabel('Ψ ');
h1=legend('1:v=16km/h','2:v=36km/h','3:v=128km/h');
set(h1,'Box','off');            %不显示方框
grid on

figure
semilogx(w,20*log10(Djw1),w,20*log10(Djw2),w,20*log10(Djw3));
ylabel('dB ');
xlabel('Ψ ');
h1=legend('1:v=16km/h','2:v=36km/h','3:v=128km/h');
set(h1,'Box','off');            %不显示方框
grid on

figure
semilogx(w,20*log10(Djw1.*Fs1),w,20*log10(Djw2.*Fs2),w,20*log10(Djw3.*Fs3));
ylabel('dB ');
xlabel('Ψ ');
h1=legend('1:v=16km/h','2:v=36km/h','3:v=128km/h');
set(h1,'Box','off');            %不显示方框
grid on

%% 角频率是啥意思？




