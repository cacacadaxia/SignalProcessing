
close all;
clear all;
fs = 10;
t = 0:1/fs:1;
f = 1;
x = sin(2*pi*f*t);


N = 1000;
f = (-N/2:N/2-1)/N*fs;
x_fft = fftshift(fft(x,N));
figure;
plot(f , abs(x_fft) )