%¶þ½×²¹³¥ÂË²¨Æ÷ ¼¶Áª ¿¹»ìÂË²¨Æ÷

%¿¹»ìÂË²¨Æ÷
k = 10^5 / 2^14;
a = k/2;
b = (k^2-a^2)^(1/2);

v   = 16;
v   = v/3.6;
pesi= 0.0001 : 0.0001 : 10;
Fs1 = (a^2+b^2) ./(-(2*pi*pesi*v).^2 + 2*1j*a*2*pi.*pesi*v + (a^2+b^2));
v   = 72;
v   = v/3.6;
Fs2 = (a^2+b^2) ./(-(2*pi*pesi*v).^2 + 2*1j*a*2*pi.*pesi*v + (a^2+b^2));
v   = 128;
v   = v/3.6;
Fs3 = (a^2+b^2) ./(-(2*pi*pesi*v).^2 + 2*1j*a*2*pi.*pesi*v + (a^2+b^2));

%% ²¹³¥ÂË²¨Æ÷
%% ÂË²¨Æ÷ÏµÊý
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
xlabel('¦· ');
h1=legend('1:v=16km/h','2:v=36km/h','3:v=128km/h');
set(h1,'Box','off');            %²»ÏÔÊ¾·½¿ò
grid on

figure
semilogx(w,20*log10(Djw1),w,20*log10(Djw2),w,20*log10(Djw3));
ylabel('dB ');
xlabel('¦· ');
h1=legend('1:v=16km/h','2:v=36km/h','3:v=128km/h');
set(h1,'Box','off');            %²»ÏÔÊ¾·½¿ò
grid on


%% »¥²¹ÂË²¨Æ÷¼ÓÉÏ¿¹»ìµþÂË²¨Æ÷
figure
semilogx(w,20*log10(Djw1.*Fs1),w,20*log10(Djw2.*Fs2),w,20*log10(Djw3.*Fs3));
ylabel('dB ');
xlabel('¦· ');
h1=legend('1:v=16km/h','2:v=36km/h','3:v=128km/h');
set(h1,'Box','off');            %²»ÏÔÊ¾·½¿ò
grid on





