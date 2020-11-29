
% =========================================================================
%
%                  ������˲����ķ���
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 11��9��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.һ�׵�ͨ�˲�����ȥ�Ʊ��˲�����Ϊ�о��Ķ���
%        2.
%       3. 
%--------------------------------------------------------------------------

%% ģ���˲����������˲�����ת��
%%����ͨ�õ������źŴ����ģ��
clear;
close all;
%һ��ģ�⿹���˲���
omiga1 = 10^5 / 2^17;
f  = 0.001 : 0.0001 :100;   % ģ��Ƶ��Hz
omiga =  2*pi.*f;           % ģ���Ƶ�� rad/sec
Bs = omiga1 ./ (1j.*(omiga) + omiga1);
figure1 = figure('Color',[1 1 1]);
semilogx( f , 20*log10(abs(Bs)) ,'LineWidth',1);
xlabel('ģ��Ƶ�� (Hz)');
ylabel('dB');            %��ֹƵ��Ϊ 0.76/2/pi = 0.121Hz
title('һ�׿����˲���');
grid on
f = 0.01:0.01:250;%%����500�ͺܹ֡�
fs = 500;	% ����Ĳ���Ƶ��Ϊ500Hz

z_1 = exp( - 1j*2*pi*f/fs);%%f����
% ��������
Omega_1 = (10^5)/(2^17);
v = 100/3.6;
T = 1/fs;
b = [ Omega_1 * T 0];
a = [ 1 + Omega_1*T ,  (-1)];
% ��ͼ
Bz = Omega_1 * T./ ( 1 + Omega_1*T - z_1 );
hold on
semilogx( f , 20*log10(abs(Bz)),'LineWidth',1);
legend ģ���˲��� �����˲���
set(gca,'Fontname','Times New Roman','fontsize',14); 
%% ��֤���ۣ�
% % �����ٶȶ��٣����չ̶�ʱ������������˲�����ģ���˲����ڿռ����ϵı���һ��
% ���仰˵���ܲ����ʶ��٣������˲�����Ӧ�ñ���һ��
% clear;
% % close all;
% %һ��ģ�⿹���˲���
% omiga1 = 10^5 / 2^17;
% f  = 0.001 : 0.0001 :100;   % ģ��Ƶ��Hz
% omiga =  2*pi.*f;           % ģ���Ƶ�� rad/sec
% Bs = omiga1 ./ (1j.*(omiga) + omiga1);
% figure1 = figure('Color',[1 1 1]);
% v = 100/3.6;
% semilogx( v./f , 20*log10(abs(Bs)) ,'LineWidth',1);
% xlabel('���� /m');
% ylabel('dB');            %��ֹƵ��Ϊ 0.76/2/pi = 0.121Hz
% title('һ�׿����˲���');
% grid on
% f = 0.01:0.01:250;%%����500�ͺܹ֡�
% fs = 500;	% ����Ĳ���Ƶ��Ϊ500Hz
% z_1 = exp( - 1j*2*pi*f/fs);%%f����
% 
% % ��������
% Omega_1 = (10^5)/(2^17);
% T = 1/fs;
% b = [ Omega_1 * T 0];
% a = [ 1 + Omega_1*T ,  (-1)];
% % ��ͼ
% Bz = Omega_1 * T./ ( 1 + Omega_1*T - z_1 );
% hold on
% semilogx(v./f , 20*log10(abs(Bz)),'LineWidth',1);
% legend ģ���˲��� �����˲���
% set(gca,'Fontname','Times New Roman','fontsize',14);
%% �����˲�����ȥ�Ʊ��˲�����������
clear all;
close all;
f = 0.01:0.01:250;%%����500�ͺܹ֡�
fs = 500;	% ����Ĳ���Ƶ��Ϊ500Hz
% ��������
Omega_1 = (10^5)/(2^17);%%��ֹƵ��
% Omega_1 = 0.9;%%��ֹƵ��
omiga =  2*pi.*f;           % ģ���Ƶ�� rad/sec
Bs = Omega_1 ./ (1j.*(omiga) + Omega_1);

%% 50km/h
figure1 = figure('Color',[1 1 1]);
v1 = 50/3.6;
semilogx(v1./f , 20*log10(abs(Bs)),'LineWidth',1);hold on;
Ti = 0.25/v1;
z_1 = exp( - 1j*2*pi* f/v1 *0.25);
Cz1 = ( exp(Omega_1*Ti/2) - exp(- Omega_1*Ti/2)*z_1 )/(Omega_1*Ti);
hold on;semilogx(v1./f , 20*log10(abs(Cz1)),'LineWidth',1);

%% 120km/h
v2 = 120/3.6;
semilogx(v2./f , 20*log10(abs(Bs)));
Ti = 0.25/v2;
z_1 = exp( - 1j*2*pi* f/v2 *0.25);
Cz2 = ( exp(Omega_1*Ti/2) - exp(- Omega_1*Ti/2)*z_1 )/(Omega_1*Ti);
hold on;semilogx(v2./f , 20*log10(abs(Cz2)),'LineWidth',1);

%% 350km/h
v3 = 200/3.6;
semilogx(v3./f , 20*log10(abs(Bs)),'LineWidth',1);
T = 0.25/v3;
z_1 = exp( - 1j*2*pi* f/v3 * 0.25);
Cz3 = ( exp(Omega_1*T/2) - exp(- Omega_1*T/2)*z_1 )/(Omega_1*T);
hold on;semilogx(v3./f , 20*log10(abs(Cz3)),'LineWidth',1);
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('\lambda /m')
ylabel('Mag /dB')
legend 50km/h 50km/h 120km/h 120km/h 350km/h 350km/h
% legend 50km/h 120km/h 350km/h
grid on;

%% Cz��freqz�����ǲ���һ����
% figure;
% semilogx(v3./f , 20*log10(abs(Cz)) , 'r','LineWidth',1);
% b = [1+Omega_1*T/2 , -1+Omega_1*T/2];
% a =  Omega_1*T;
% [Czz ,f] = freqz(b,a , length(f),1/T);
% hold on;
% semilogx(v3./f , 20*log10(abs(Czz)),'k','LineWidth',1);

%%
figure1 = figure('Color',[1 1 1]);
semilogx(v1./f , 20*log10(abs(Bs.*Cz1)),'LineWidth',1);
hold on;
semilogx(v2./f , 20*log10(abs(Bs.*Cz2)),'LineWidth',1);
semilogx(v3./f , 20*log10(abs(Bs.*Cz3)),'LineWidth',1);
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('\lambda /m')
ylabel('Mag /dB')
title('������˲���+ȥ�Ʊ��˲���')
legend 50km/h 120km/h 200km/h

