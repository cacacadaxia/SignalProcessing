close all;
clear all;
%һ��ģ�⿹���˲���
omiga1 = 10^5 / 2^17;
f  = 0.001 : 0.0001 :100;  % ģ��Ƶ��Hz
omiga =  2*pi.*f;         % ģ���Ƶ�� rad/sec
H = omiga1 ./ (1j.*(omiga) + omiga1);
figure
semilogx( f , 20*log10(H) );
xlabel('ģ��Ƶ�� (Hz)');
ylabel('dB');            % ��ֹƵ��Ϊ 0.76/2/pi = 0.121Hz
title('һ�׿����˲���');


%%
%��ͬ�ٶ��£�һ��ģ�⿹���˲�����Ƶ����
%�ٶ�Ϊ16km/h 36km/h 128km/h
lamda = 1:0.001:10000;        %����
pesi  = 1./lamda;             %�ռ�Ƶ��
v  = 16; 
v = v/3.6;
H1 = omiga1 ./ (1j*(2*pi*v.*pesi) + omiga1);
v  = 36;
v = v/3.6;
H2 = omiga1 ./ (1j*(2*pi*v.*pesi) + omiga1);
v  = 128;
v = v/3.6;
H3 = omiga1 ./ (1j*(2*pi*v.*pesi) + omiga1);

figure
semilogx(lamda,20*log10(H1),lamda,20*log10(H2),lamda,20*log10(H3));
semilogx(pesi,20*log10(H1),pesi,20*log10(H2),pesi,20*log10(H3));
title ('�����˲����ռ���Ƶ������')
xlabel('�ռ�Ƶ��(1/m)' );
ylabel('��ֵ(dB)');
legend('16km/h','36km/h','128km/h');

%% �����˲���




