
% =========================================================================
%
%                  抗混叠滤波器的方法
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 11月9日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.一阶低通滤波器和去移变滤波器作为研究的对象
%        2.
%       3. 
%--------------------------------------------------------------------------

%% 模拟滤波器与数字滤波器的转换
%%这是通用的数字信号处理的模型
clear;
close all;
%一阶模拟抗混滤波器
omiga1 = 10^5 / 2^17;
f  = 0.001 : 0.0001 :100;   % 模拟频率Hz
omiga =  2*pi.*f;           % 模拟角频率 rad/sec
Bs = omiga1 ./ (1j.*(omiga) + omiga1);
figure1 = figure('Color',[1 1 1]);
semilogx( f , 20*log10(abs(Bs)) ,'LineWidth',1);
xlabel('模拟频率 (Hz)');
ylabel('dB');            %截止频率为 0.76/2/pi = 0.121Hz
title('一阶抗混滤波器');
grid on
f = 0.01:0.01:250;%%到了500就很怪。
fs = 500;	% 惯组的采样频率为500Hz

z_1 = exp( - 1j*2*pi*f/fs);%%f代表
% 参数配置
Omega_1 = (10^5)/(2^17);
v = 100/3.6;
T = 1/fs;
b = [ Omega_1 * T 0];
a = [ 1 + Omega_1*T ,  (-1)];
% 画图
Bz = Omega_1 * T./ ( 1 + Omega_1*T - z_1 );
hold on
semilogx( f , 20*log10(abs(Bz)),'LineWidth',1);
legend 模拟滤波器 数字滤波器
set(gca,'Fontname','Times New Roman','fontsize',14); 
%% 验证结论：
% % 不论速度多少，按照固定时间采样的数字滤波器域模拟滤波器在空间域上的表现一致
% 换句话说不管采样率多少，两个滤波器都应该表现一致
% clear;
% % close all;
% %一阶模拟抗混滤波器
% omiga1 = 10^5 / 2^17;
% f  = 0.001 : 0.0001 :100;   % 模拟频率Hz
% omiga =  2*pi.*f;           % 模拟角频率 rad/sec
% Bs = omiga1 ./ (1j.*(omiga) + omiga1);
% figure1 = figure('Color',[1 1 1]);
% v = 100/3.6;
% semilogx( v./f , 20*log10(abs(Bs)) ,'LineWidth',1);
% xlabel('波长 /m');
% ylabel('dB');            %截止频率为 0.76/2/pi = 0.121Hz
% title('一阶抗混滤波器');
% grid on
% f = 0.01:0.01:250;%%到了500就很怪。
% fs = 500;	% 惯组的采样频率为500Hz
% z_1 = exp( - 1j*2*pi*f/fs);%%f代表
% 
% % 参数配置
% Omega_1 = (10^5)/(2^17);
% T = 1/fs;
% b = [ Omega_1 * T 0];
% a = [ 1 + Omega_1*T ,  (-1)];
% % 画图
% Bz = Omega_1 * T./ ( 1 + Omega_1*T - z_1 );
% hold on
% semilogx(v./f , 20*log10(abs(Bz)),'LineWidth',1);
% legend 模拟滤波器 数字滤波器
% set(gca,'Fontname','Times New Roman','fontsize',14);
%% 数字滤波器和去移变滤波器参数配置
clear all;
close all;
f = 0.01:0.01:250;%%到了500就很怪。
fs = 500;	% 惯组的采样频率为500Hz
% 参数配置
Omega_1 = (10^5)/(2^17);%%截止频率
% Omega_1 = 0.9;%%截止频率
omiga =  2*pi.*f;           % 模拟角频率 rad/sec
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

%% Cz用freqz来画是不是一样的
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
title('抗混叠滤波器+去移变滤波器')
legend 50km/h 120km/h 200km/h

