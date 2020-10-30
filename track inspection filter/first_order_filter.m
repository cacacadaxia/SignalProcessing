close all
clear 

%һ��ģ�⿹���˲���
omiga1 = 10^5 / 2^17;
f  = 0.001 : 0.0001 :100;  % ģ��Ƶ��Hz
omiga =  2*pi.*f;         % ģ���Ƶ�� rad/sec
H = omiga1 ./ (1j.*(omiga) + omiga1);
figure1 = figure('Color',[1 1 1]);
semilogx( f , 20*log10(H) );
xlabel('ģ��Ƶ�� (Hz)');
ylabel('dB');            % ��ֹƵ��Ϊ 0.76/2/pi = 0.121Hz
title('һ�׿����˲���');


%% �����ٶȵ�Ӱ��
%��ͬ�ٶ��£�һ��ģ�⿹���˲�����Ƶ����
%�ٶ�Ϊ16km/h 36km/h 128km/h
lamda = 1:0.001:10000;        %����
pesi  = 1./lamda;             %�ռ�Ƶ��
v  = 16; v = v/3.6;
H1 = omiga1 ./ (1j*(2*pi*v.*pesi) + omiga1);
v  = 36;v = v/3.6;
H2 = omiga1 ./ (1j*(2*pi*v.*pesi) + omiga1);
v  = 128;v = v/3.6;
H3 = omiga1 ./ (1j*(2*pi*v.*pesi) + omiga1);

figure1 = figure('Color',[1 1 1]);
semilogx(lamda,20*log10(H1),lamda,20*log10(H2),lamda,20*log10(H3));
% semilogx(pesi,20*log10(H1),pesi,20*log10(H2),pesi,20*log10(H3));
title ('�����˲����ռ���Ƶ������')
xlabel('�ռ�Ƶ��(1/m)' );
ylabel('��ֵ(dB)');
legend('16km/h','36km/h','128km/h');
 
%% ����ķ���ֵ��ѧϰһ��
%% s = 1j*(2*pi*v.*pesi)
%% z^-1 = exp(-1j*2*pi.*pesi*v*Ti)
%% z = exp(1j*2*pi.*pesi*v*Ti);
%% 
%% ���������ٶ�v�ı��ʽ
%% ���û��v�������£�
%% �ο����ģ������˲������ڹ������е�Ӧ��
%% 
% % һ�������˲���
% % ��ģ���˲���ת��Ϊ�����˲��� �����ַ�
% T = 0.002;
% v = 16;
% H4 = omiga1*T  ./ (1+omiga1*T-exp(-1j.*2*pi*v.*pesi*T));
% v  = 36;
% H5 = omiga1*T  ./ (1+omiga1*T-exp(-1j.*2*pi*v.*pesi*T));
% v  = 128;
% H6 = omiga1*T  ./ (1+omiga1*T-exp(-1j.*2*pi*v.*pesi*T));
% figure
% semilogx(lamda,20*log10(H4),lamda,20*log10(H5),lamda,20*log10(H6));
% title ('���ֿ����˲���');
% xlabel('����(m)' );
% ylabel('��ֵ(dB)');
% legend('16km/h','36km/h','128km/h');


%%
%һ�����ֲ����˲���
% ��Խ���̵�ϵ�������ڼ�
v    = 16;
v    = v/3.6;
Ti   = 0.25/v;
par1 = omiga1 * Ti;
B1 = ( (1+par1/2)-(1-par1/2).*exp(-1j*2*pi.*pesi*v*Ti) ) ./ par1;

v    = 36;
v    = v/3.6;
Ti   = 0.25/v;
par1 = omiga1 * Ti;
B2 = ( (1+par1/2)-(1-par1/2).*exp(-1j*2*pi.*pesi*v*Ti) ) ./ par1;

v    = 128;
v    = v/3.6;
Ti   = 0.25/v;
par1 = omiga1 * Ti;
B3 = ( (1+par1/2)-(1-par1/2).*exp(-1j*2*pi.*pesi*v*Ti) ) ./ par1;

% figure;
hold on
% semilogx(pesi,20*log10(B1),pesi,20*log10(B2),pesi,20*log10(B3));
semilogx(lamda,20*log10(B1),lamda,20*log10(B2),lamda,20*log10(B3));
title ('�����˲����ռ���Ƶ������');
xlabel('�ռ�Ƶ��(1/m)' );
ylabel('��ֵ(dB)');
legend('16km/h','36km/h','128km/h');


%%
%ģ�⿹���˲���+���ֲ����˲���
figure1 = figure('Color',[1 1 1]);
% semilogx(pesi,20*log10(B1.*H1),pesi,20*log10(B2.*H2),pesi,20*log10(B3.*H3));
semilogx(lamda,20*log10(B1.*H1),lamda,20*log10(B2.*H2),lamda,20*log10(B3.*H3));
title ('ģ�⿹�졢���ֲ����˲����ռ���Ƶ������');
% xlabel('�ռ�Ƶ��(1/m)' );
xlabel('����(m)' );
ylabel('��ֵ(dB)');
legend('16km/h','36km/h','128km/h');

%% ע��
%%ͬ��Ҳ��һ����ͨ�˲�������1m���ϵĲ�������������


% %���ֿ����˲���+���ֲ����˲���
% figure
% semilogx(lamda,20*log10(B1.*H4),lamda,20*log10(B2.*H5),lamda,20*log10(B3.*H6));
% title ('���ֿ����˲���+���ֲ����˲���');
% xlabel('����(m)' );
% ylabel('��ֵ(dB)');
% legend('16km/h','36km/h','128km/h');




