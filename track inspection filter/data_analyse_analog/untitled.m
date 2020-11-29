clear all;
Fs=1000;
t=0:1/Fs:.3;
f1 = 200;
sig=cos(2*pi*t*f1)+randn(size(t));
Hs=spectrum.welch;
psd(Hs,sig,'Fs',Fs);

%%
N = length(sig);
x = (-N/2:N/2-1)/N*Fs;
figure;plot(x,20*log10(abs(fftshift(fft(sig)))));