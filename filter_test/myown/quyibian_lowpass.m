

% =========================================================================
%
%                  抗混叠滤波器的方法
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月16日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.一阶低通滤波器和去移变滤波器作为研究的对象
%        2.
%       3. 
%--------------------------------------------------------------------------

%% 模拟滤波器与数字滤波器


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
v = 50/3.6;
T = 1/fs;
b = [ Omega_1 * T 0];
a = [ 1 + Omega_1*T ,  (-1)];
% 画图
Bz = Omega_1 * T./ ( 1 + Omega_1*T - z_1 );
hold on
semilogx(f , 20*log10(abs(Bz)),'LineWidth',1);
legend 模拟滤波器 数字滤波器
set(gca,'Fontname','Times New Roman','fontsize',14);

 
%% 不论速度多少，按照固定时间采样的数字滤波器域模拟滤波器在空间域上的表现一致
clear;
% close all;
%一阶模拟抗混滤波器
omiga1 = 10^5 / 2^17;
f  = 0.001 : 0.0001 :100;   % 模拟频率Hz
omiga =  2*pi.*f;           % 模拟角频率 rad/sec
Bs = omiga1 ./ (1j.*(omiga) + omiga1);
figure1 = figure('Color',[1 1 1]);
v = 100/3.6;
semilogx( v./f , 20*log10(abs(Bs)) ,'LineWidth',1);
xlabel('波长 /m');
ylabel('dB');            %截止频率为 0.76/2/pi = 0.121Hz
title('一阶抗混滤波器');
grid on
f = 0.01:0.01:250;%%到了500就很怪。
fs = 500;	% 惯组的采样频率为500Hz
z_1 = exp( - 1j*2*pi*f/fs);%%f代表

% 参数配置
Omega_1 = (10^5)/(2^17);
T = 1/fs;
b = [ Omega_1 * T 0];
a = [ 1 + Omega_1*T ,  (-1)];
% 画图
Bz = Omega_1 * T./ ( 1 + Omega_1*T - z_1 );
hold on
semilogx(v./f , 20*log10(abs(Bz)),'LineWidth',1);
legend 模拟滤波器 数字滤波器
set(gca,'Fontname','Times New Roman','fontsize',14);


%% 数字滤波器和去移变滤波器参数配置
clear all;
f = 0.01:0.01:250;%%到了500就很怪。
fs = 500;	% 惯组的采样频率为500Hz
z_1 = exp( - 1j*2*pi*f/fs);%%f代表
% 参数配置
Omega_1 = (10^5)/(2^17);%%截止频率
figure1 = figure('Color',[1 1 1]);

%% 50km/h
v = 50/3.6;
T = 1/fs;
b = [ Omega_1 * T , 0];
a = [ 1 + Omega_1*T ,  (-1)];
% 画图
z_1 = exp( - 1j*2*pi*f/fs);%%f代表
Bz = Omega_1 * T./ ( 1 + Omega_1*T - z_1 );
semilogx(v./f , 20*log10(abs(Bz)),'LineWidth',1);
T = 0.25/v;
z_1 = exp( - 1j*2*pi* f/v *0.25);
Cz = ( exp(Omega_1*T/2) - exp(- Omega_1*T/2)*z_1 )/(Omega_1*T);
hold on;semilogx(v./f , 20*log10(abs(Cz)));
Bz1 = Bz;
Cz1 = Cz;


%% 120km/h
v = 120/3.6;
T = 1/fs;
b = [ Omega_1 * T , 0];
a = [ 1 + Omega_1*T ,  (-1)];
% 画图
z_1 = exp( - 1j*2*pi*f/fs);%%f代表
Bz = Omega_1*T./ ( 1 + Omega_1*T - z_1 );
hold on
semilogx(v./f , 20*log10(abs(Bz)));
T = 0.25/v;
z_1 = exp( - 1j*2*pi* f/v *0.25);
Cz = ( exp(Omega_1*T/2) - exp(- Omega_1*T/2)*z_1 )/(Omega_1*T);
hold on;semilogx(v./f , 20*log10(abs(Cz)));

Bz2 = Bz;
Cz2 = Cz;

%% 350km/h
v = 350/3.6;
T = 1/fs;
b = [ Omega_1 * T , 0];
a = [ 1 + Omega_1*T ,  (-1)];
% 画图
z_1 = exp( - 1j*2*pi*f/fs);%%f代表
Bz = Omega_1*T./ ( 1 + Omega_1*T - z_1 );
semilogx(v./f , 20*log10(abs(Bz)));
T = 0.25/v;
z_1 = exp( - 1j*2*pi* f/v * 0.25);
Cz = ( exp(Omega_1*T/2) - exp(- Omega_1*T/2)*z_1 )/(Omega_1*T);
hold on;semilogx(v./f , 20*log10(abs(Cz)));
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('\lambda /m')
ylabel('Mag /dB')
legend 50km/h 50km/h 120km/h 120km/h 350km/h 350km/h
% legend 50km/h 120km/h 350km/h
grid on;

Bz3 = Bz;
Cz3 = Cz;

%%
figure1 = figure('Color',[1 1 1]);
semilogx(50/3.6./f , 20*log10(abs(Bz1.*Cz1)));
hold on;
semilogx(120/3.6./f , 20*log10(abs(Bz2.*Cz2)));
semilogx(350/3.6./f , 20*log10(abs(Bz3.*Cz3)));
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('\lambda /m')
ylabel('Mag /dB')
title('抗混叠滤波器+去移变滤波器')
legend 50km/h 120km/h 350km/h


%% 模拟滤波器+去移变滤波器
