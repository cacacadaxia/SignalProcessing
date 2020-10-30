

close all;
omega = 0.76;
figure;freqs(omega,[1,omega]);

[h1,f1] = freqs(omega,[1,omega]);
figure;semilogx(f1/2/pi,(h1));
fs   = 1000;
[h2,f2] = freqz([omega/fs],[1+omega/fs , -1],10000,fs);
hold on;semilogx(f2,h2);




