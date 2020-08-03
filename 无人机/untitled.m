close all;
% 普通的低通滤波器
clear all;
IMU_SAMPLE_RATE = 100;
FILTER_CUTOFF = 30;

sample = IMU_SAMPLE_RATE;
cutoff = FILTER_CUTOFF;

% function SetFilter(sample,cutoff)
fr = sample/cutoff;
ohm = tan(pi/fr);
c = 1+2*cos(pi/4)*ohm + ohm*ohm;

if cutoff>0
   b01 = ohm*ohm/c ;
   b11 = 2*b01;
   b21 = b01;
   a11 = 2*(ohm*ohm-1)/c;
   a21 = (1-2*cos(pi/4)*ohm + ohm^2)/c;
end
% end

b = [b01,b11,b21];
a = [0 , a11,a21];
[ H , f ] = freqz(b,a,512,sample);%%f为sample的一半
figure;
plot(f, 20*log10(abs(H)),'LineWidth',1);
