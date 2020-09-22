a = [1 0.4 1];

b = [0.2 0.3 1];

w = logspace(-1, 1);

freqs(b, a, w);

h = freqs(b, a, w);
mag = abs(h);
phase = angle(h);
figure;
subplot(2,1,1), loglog(w,mag)
subplot(2,1,2), semilogx(w,phase);
grid on

