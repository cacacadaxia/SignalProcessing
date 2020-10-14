close all
clear
clc


%% �ǲ��Ǿ���Fs��Gz�����˲�����
%% 
% ����butterworth�����˲���
% ʱ���Ƶ����
a = 10^5 / 2^14;
f = 0.0001 : 0.0001 :10;        % ģ��Ƶ��Hz
omiga = 2 * pi .* f;              % ģ���Ƶ�� rad/sec
Fs1   = a^2 ./((1j.*omiga).^2 + 1j.*omiga*a + a^2);
figure;
semilogx(f,20*log10(Fs1));
xlabel('��(Hz)');
ylabel('dB ');
title('���׿����˲���ʱ���Ƶ��������');

% ����butterworth�����˲���
% �ռ����Ƶ����
k = 10^5 / 2^14;
a = k/2;
b = (k^2-a^2)^(1/2);

%�ռ�Ƶ��
pesi = 0.0001 : 0.0001 : 10;     %�ռ䲨�� 1m ~ 1000m
v   = 16/3.6;
Fs1 = (a^2+b^2) ./(-(2*pi*pesi*v).^2 + 2*1j*a*2*pi.*pesi*v + (a^2+b^2));
v   = 36/3.6;
Fs2 = (a^2+b^2) ./(-(2*pi*pesi*v).^2 + 2*1j*a*2*pi.*pesi*v + (a^2+b^2));
v   = 128/3.6;
Fs3 = (a^2+b^2) ./(-(2*pi*pesi*v).^2 + 2*1j*a*2*pi.*pesi*v + (a^2+b^2));


figure
semilogx(pesi,20*log10(Fs1),pesi,20*log10(Fs2),pesi,20*log10(Fs3));
ylabel('dB ');
xlabel('��(1/m) ');
title('���׿����˲����ռ����Ƶ��������');
h1=legend('1:v=16km/h','2:v=36km/h','3:v=128km/h');
set(h1,'Box','off');            



%% ���׿����˲�����Ӧ�Ĳ����˲���
%% Ӧ��������������˲���
% �����˲���
w = 0.0001 : 0.0001 :10;           %�ռ䲨�� 1m ~ 1000m
v = 16/3.6;
t = 0.25/v;
a1 = ( 1 + a*t + a^2*t*t )/( t^2*(a^2+b^2) );
a2 = 2*(1-a^2*t^2) / ( t^2*(a^2+b^2) ) ;
a3 = ( 1 - a*t + a^2*t*t ) / ( t^2*(a^2+b^2) ) ;
Djw1 = a1 - a2*exp(-1j*2*pi*0.25.*w) + a3*exp(-2j*2*pi*0.25.*w);

v  = 36/3.6;
t  = 0.25/v;
a1 = ( 1 + a*t + a^2*t*t )/( t^2*(a^2+b^2) );
a2 = 2*(1-a^2*t^2) / ( t^2*(a^2+b^2) ) ;
a3 = ( 1 - a*t + a^2*t*t ) / ( t^2*(a^2+b^2) ) ;
Djw2 = a1-a2*exp(-1j*2*pi*0.25.*w)+a3*exp(-2j*2*pi*0.25.*w);

v  = 128/3.6;
t  = 0.25/v;
a1 = ( 1 + a*t + a^2*t*t )/( t^2*(a^2+b^2) );
a2 = 2*(1-a^2*t^2) / ( t^2*(a^2+b^2) ) ;
a3 = ( 1 - a*t + a^2*t*t ) / ( t^2*(a^2+b^2) ) ;
Djw3 = a1-a2*exp(-1j*2*pi*0.25.*w)+a3*exp(-2j*2*pi*0.25.*w);

figure
semilogx(w,20*log10(Djw1),w,20*log10(Djw2),w,20*log10(Djw3));
xlabel('����1/m�� ');
ylabel('dB ');
title('���ײ����˲����ռ����Ƶ��������');
h1=legend('1:v=16km/h','2:v=36km/h','3:v=128km/h');
set(h1,'Box','off');            

figure;
semilogx(w,20*log10(Djw1.*Fs1),w,20*log10(Djw2.*Fs2),'b',w,20*log10(Djw3.*Fs3));
xlabel('����1/m�� ');
ylabel('dB ');
title('����-�����˲���������Ƶ��������');
h1=legend('1:v=16km/h','2:v=36km/h','3:v=128km/h');
set(h1,'Box','off');           %����ʾ����




%% ����˲��������˿����+�����˲�����Ч������϶��ɵĵ�Ƶ�˲���������Ч�Ŀ����
%%�����²���Ϊ1m���ϵģ����ⲿ������������Ҫ��
%%����ԱȲ�һ���ԣ�������ʱ��ȥ���Ǿ�����

figure;semilogx(1./w,20*log10(Djw1.*Fs1),1./w,20*log10(Djw2.*Fs2),'b',1./w,20*log10(Djw3.*Fs3));
xlabel('\lamda \m');
%%���Ի��ǵ�Ƶ�˲���




%%
% hold on;
% a = 10^5 / 2^14;
% f = 0.0001 : 0.0001 :10;        % ģ��Ƶ��Hz
% omiga = 2 * pi .* f;              % ģ���Ƶ�� rad/sec
% Fs1   = a^2 ./((1j.*omiga).^2 + 1j.*omiga*a + a^2);
% semilogx(1./f,20*log10(Fs1));
% legend;
% xlabel('WaveLength m')
% ylabel('Mag dB')


% %%
% figure;
% %�ռ�Ƶ��
% pesi = 0.0001 : 0.0001 : 10;     %�ռ䲨�� 1m ~ 1000m
% v   = 3.6/3.6;
% Fs1 = (a^2+b^2) ./(-(2*pi*pesi*v).^2 + 2*1j*a*2*pi.*pesi*v + (a^2+b^2));
% semilogx(1./w,20*log10(Fs1));
% hold on;
% a = 10^5 / 2^14;
% f = 0.0001 : 0.0001 :10;        % ģ��Ƶ��Hz
% omiga = 2 * pi .* f;              % ģ���Ƶ�� rad/sec
% Fs1   = a^2 ./((1j.*omiga).^2 + 1j.*omiga*a + a^2);
% semilogx(1./f,20*log10(Fs1));
% legend;
% xlabel('WaveLength m')
% ylabel('Mag dB')
