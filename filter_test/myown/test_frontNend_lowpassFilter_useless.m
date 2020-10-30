


close all;
clear all;
omega_0 = 5.71;
omega_0 = 6.67;
RC = 1/omega_0;

K = 2;
b = K;
a = [RC^2 , 3*RC-K*RC,1];
figure;
[H,tmp ] = freqs(b,a);
semilogx(tmp/pi/2 , 20*log10(H),'LineWidth',2);
grid on;


%%
% ----- ¶þ½×¿¹»ìµþ + ²¹³¥ -------------
clear all;
clc;
% ¶þ½×µÍÍ¨ÂË²¨Æ÷
w2 = (10^5)/(2^14);
v1=  440/3.6;
t1= 0.25/v1;

w2t1 = w2*t1;
b2 = [(w2t1)^2];
a2 = [1+w2t1+(w2t1)^2  ,- (2 + w2t1)  ,1];
[h2 f2] = freqz(b2,a2,800000,500);
figure;suptitle ('¶þ½×Êý×Ö¿¹»ìµþÂË²¨Æ÷ºÍ²¹³¥ÂË²¨Æ÷');
semilogx(v1./f2,20*log10(abs(h2)));hold on;

%%
w2 = (10^5)/(2^14);
b = [(w2)^2];
a = [1,w2,w2^2];
figure;
[H,f] = freqs(b,a);
semilogx(f/2/pi,20*log10(abs(H)));hold on;
semilogx(f2,20*log10(abs(h2)));hold on;
legend 1 2
