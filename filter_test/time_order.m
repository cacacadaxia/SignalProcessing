%һ�׶����˲���ʱ���ϵĸ�Ƶ����
% ��Ϊ�趨ʱ��Ƶ��Ϊ  t= 0.002
clear all;
close all;
% һ��ģ��
Q1 = 2^17;
Q1 = 100000/Q1;
b = [0 Q1];
a =[1 Q1];
[h0 f0] = freqs(b,a);
figure();suptitle('һ��ģ�⿹���ʱ��Ƶ��������');subplot(2,1,1);
semilogx(f0/(2*pi),20*log10(abs(h0)));
xlabel('Ƶ�ʣ�Hz��');ylabel('��ֵ(dB)');
subplot(2,1,2);
semilogx(f0/(2*pi),angle(h0)*180/pi);
xlabel('Ƶ�ʣ�Hz��');ylabel('��λ');
% һ��ʱ����
Q1 = 2^17;
Q1 = 100000/Q1;
t = 0.002;
b = [Q1*t 0];
a = [1+Q1*t -1];
[h,f] = freqz (b,a,15000,500);
figure();suptitle('һ�����ֿ����ʱ��Ƶ��������');subplot(2,1,1);
semilogx(f,20*log10(abs(h)));
xlabel('Ƶ�ʣ�Hz��');ylabel('��ֵ(dB)');
subplot(2,1,2);
semilogx(f,angle(h)*180/pi);
xlabel('Ƶ�ʣ�Hz��');ylabel('��λ');
% ����
bc = [1+(Q1*t/2)  (-1)+(Q1*t/2)];
ac = [Q1*t 0];
bc = [ 1+(Q1*t/2)  (-1)+(Q1*t/2)];
ac = [Q1*t 0];
[h1,f1] = freqz (bc,ac,15000,500);
figure();suptitle('һ�����ֲ���ʱ��Ƶ��������');subplot(2,1,1);
semilogx(f1,20*log10(abs(h1)));
xlabel('Ƶ�ʣ�Hz��');ylabel('��ֵ(dB)');
subplot(2,1,2);
semilogx(f1,angle(h1)*180/pi);
xlabel('Ƶ�ʣ�Hz��');ylabel('��λ');

% ----����ģ��-----
Q = 2^14;
Q = 100000/Q;
b = [0 0 Q^2];
a = [1 Q Q^2];
[h0 f0] = freqs(b,a);
figure();suptitle('����ģ�⿹���ʱ��Ƶ��������');subplot(2,1,1);
semilogx(f0/(2*pi),20*log10(abs(h0)));
xlabel('Ƶ�ʣ�Hz��');ylabel('��ֵ(dB)');
subplot(2,1,2);
semilogx(f0/(2*pi),angle(h0)*180/pi);
xlabel('Ƶ�ʣ�Hz��');ylabel('��λ');
% ------����ʱ���ϵ�-----
Q = 2^14;
Q = 100000/Q;
t = 0.002;
b0 = (Q*t)^2;
a0 = (Q*t)^2+Q*t+1;
a1 = -(Q*t+2);
a2 = 1;
%b = [0 0 b0];
%a = [a2 a1 a0];
b = [b0 0 0];
a = [a0 a1 a2];
[h,f] = freqz (b,a,1500,500);
figure();suptitle('�������ֿ����ʱ��Ƶ��������');subplot(2,1,1);
semilogx(f,20*log10(abs(h)),'r');
xlabel('Ƶ�ʣ�Hz��');ylabel('��ֵ(dB)');
subplot(2,1,2);
semilogx(f,angle(h)*180/pi,'r');
xlabel('Ƶ�ʣ�Hz��');ylabel('��λ');
% -----����-----
a = Q/2; aa = a^2;bb = 3*(Q^2)/4;bbb = Q*sqrt(3/4);
b2 = [(1+a*t+aa*t*t)*bbb -2*(1-aa*t*t)*bbb (1-a*t+aa*t*t)*bbb];
a2 = [(aa+bb)*t*2 0 0];

b2 = [1+a*t+aa*t*t -2*(1-aa*t*t) 1-a*t+aa*t*t];
a2 = [(aa+bb)*t*t 0 0];
[h2 f2] = freqz(b2,a2,1000,500);
figure();suptitle('�������ֲ���ʱ��Ƶ��������');subplot(2,1,1);
semilogx(f2,20*log10(abs(h2)),'r');
xlabel('Ƶ�ʣ�Hz��');ylabel('��ֵ(dB)');
subplot(2,1,2);
semilogx(f2,angle(h2)*180/pi,'r');
xlabel('Ƶ�ʣ�Hz��');ylabel('��λ');

