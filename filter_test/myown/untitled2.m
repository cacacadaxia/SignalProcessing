% ½âÆ«ÂË²¨Æ÷
% ¸ßÍ¨ÂË²¨Æ÷£¿
wd = 0.001;
b = [1-wd,-1+wd];
a = [1,-1+wd];

v = 100/3.6;
t1 = 0.25/v;
[h,f] = freqz(b,a,10000,1/t1);
lamda = v./f;
figure;
semilogx(lamda,20*log10(abs(h)));


v = 200/3.6;
t1 = 0.25/v;
[h,f] = freqz(b,a,10000,1/t1);
lamda = v./f;
hold on;
semilogx(lamda,20*log10(abs(h)));
legend 1 2

