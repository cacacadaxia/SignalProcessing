%% Ò»½×Êý×Ö¿¹»ìµþÂË²¨Æ÷+²¹³¥ÂË²¨Æ÷
clear all;
clc;
%----------Ò»½×Êý×Ö¿¹»ìµþÂË²¨Æ÷------------
W1 = (10^5)/(2^17);
v1= 150/36;v2 = 200/9;v3 = 3500/36;
T1 = 0.25/v1;T2 = 0.25/v2;T3 = 0.25/v3;
%figure();suptitle('Ò»½×Êý×Ö¿¹»ìµþ')
%semilogx(f1,20*log10(abs(h1)));xlabel('ÆµÂÊ£¨Hz£©');ylabel('·ùÖµ(dB)');
figure();suptitle('Ò»½×Êý×Ö¿¹»ìµþ+²¹³¥');
b = [ W1*T1 0];
a = [ 1+W1*T1 (-1)];
[h1 f1] = freqz(b,a,800000,500);
semilogx(v1./f1,20*log10(abs(h1)));hold on;
b = [ W1*T2 0];
a = [ 1+W1*T2 (-1)];
[h1 f1] = freqz(b,a,800000,500);
semilogx(v2./f1,20*log10(abs(h1)),'g');hold on;
b = [ W1*T3 0];
a = [ 1+W1*T3 (-1)];
[h1 f1] = freqz(b,a,800000,500);
semilogx(v3./f1,20*log10(abs(h1)),'r');hold on;
xlabel('²¨³¤£¨m£©');ylabel('·ùÖµ(dB)');


%-----------Ò»½×Êý×Ö²¹³¥ÂË²¨Æ÷----------------
T1 = 0.25/v1;
bc = [ 1+(W1*T1/2)  (-1)+(W1*T1/2)];
ac = [W1*T1 0];
[h1c f1c] = freqz(bc,ac,800000,500);
semilogx(v1./f1c,20*log10(abs(h1c)));hold on;

T2 = 0.25/v2;
bc = [ 1+(W1*T2/2)  (-1)+(W1*T2/2)];
ac = [W1*T2 0];
[h1c f1c] = freqz(bc,ac,800000,500);
semilogx(v2./f1c,20*log10(abs(h1c)),'g');hold on;

bc = [ 1+(W1*T3/2)  (-1)+(W1*T3/2)];
ac = [W1*T3 0];
[h1c f1c] = freqz(bc,ac,800000,500);
semilogx(v3./f1c,20*log10(abs(h1c)),'r');hold on;

%-----------------Ò»½×Êý×Ö¿¹»ìµþ+²¹³¥ÂË²¨Æ÷-------
figure;
B1 = [2*W1*T1+(W1^2)*T1*T1 (W1^2)*T1*T1-2*W1*T1];
A1 = [2*W1*T1+2*(W1^2)*T1*T1 (-2)*W1*T1];
[H1 F1] = freqz(B1,A1,800000,500);
semilogx(v1./F1,20*log10(abs(H1)),'b');hold on;

B1 = [2*W1*T2+(W1^2)*T2*T2 (W1^2)*T2*T2-2*W1*T2];
A1 = [2*W1*T2+2*(W1^2)*T2*T2 (-2)*W1*T2];
[H1 F1] = freqz(B1,A1,800000,500);
semilogx(v2./F1,20*log10(abs(H1)),'g');hold on;

B1 = [2*W1*T3+(W1^2)*T3*T3 (W1^2)*T3*T3-2*W1*T3];
A1 = [2*W1*T3+2*(W1^2)*T3*T3 (-2)*W1*T3];
[H1 F1] = freqz(B1,A1,800000,500);
semilogx(v3./F1,20*log10(abs(H1)),'r');hold on;

xlabel('²¨³¤£¨m£©');ylabel('·ùÖµ(dB)');



