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
%  功能： 1.二阶低通滤波器和去移变滤波器作为研究的对象
%        2.已经完成了代码的编写
%       3. 
%--------------------------------------------------------------------------

%% 数字滤波器和去移变滤波器参数配置
clear all;
close all;
f = 0.01:0.01:250;%%到了500就很怪。
fs = 500;	% 惯组的采样频率为500Hz
% 参数配置
Omega_1 = (10^5)/(2^17);%%截止频率
Omega_2 = 10^5/2^14;
omiga =  2*pi.*f;           % 模拟角频率 rad/sec
% Fs = Omega_1 ./ (1j.*(omiga) + Omega_1);%%这是一阶的
Fs = Omega_2^2./( (1j.*omiga).^2 + (1j.*omiga).*Omega_2 + Omega_2^2 );

%% 50km/h
figure1 = figure('Color',[1 1 1]);
v1 = 50/3.6;
semilogx(v1./f , 20*log10(abs(Fs)),'LineWidth',1);hold on;
Ti = 0.25/v1;
z_1 = exp( - 1j*2*pi* f/v1 *0.25);
Cz1 = ( exp(Omega_1*Ti/2) - exp(- Omega_1*Ti/2)*z_1 )/(Omega_1*Ti);
a = 1/2*Omega_2; b = sqrt(1-0.5^2)*Omega_2;
a2b2 = (10^5/2^14)^2;
w2 = (10^5)/(2^14);
a = w2/2; aa = a^2;
bb = 3*(w2^2)/4;
biir = [1+a*Ti+aa*Ti*Ti ,-2*(1-aa*Ti*Ti), 1-a*Ti+aa*Ti*Ti];
% Gz1 = ( (1-z_1).^2 + (1 - z_1.^2)*a*Ti + (a*Ti + 2*a*Ti*z_1 + a*Ti*z_1.^2) ) ./ (Ti^2*a2b2);
Gz1 = (biir(1) + biir(2).*z_1 + biir(3).*z_1.^2) ./ ((aa+bb)*Ti*Ti);
hold on;semilogx(v1./f , 20*log10(abs(Gz1)));

%% 120km/h
v2 = 120/3.6;
semilogx(v2./f , 20*log10(abs(Fs)));
Ti = 0.25/v2;
z_1 = exp( - 1j*2*pi* f/v2 *0.25);
Cz2 = ( exp(Omega_1*Ti/2) - exp(- Omega_1*Ti/2)*z_1 )/(Omega_1*Ti);
% Gz2 = ( (1-z_1).^2 + (1 - z_1.^2)*a*Ti + (a*Ti + 2*a*Ti*z_1 + a*Ti*z_1.^2) ) ./ (Ti^2*a2b2);
biir = [1+a*Ti+aa*Ti*Ti ,-2*(1-aa*Ti*Ti), 1-a*Ti+aa*Ti*Ti];
Gz2 = (biir(1) + biir(2).*z_1 + biir(3).*z_1.^2) ./ ((aa+bb)*Ti*Ti);
hold on;semilogx(v2./f , 20*log10(abs(Gz2)));

%% 350km/h
v3 = 200/3.6;
v3 = 62.5;
semilogx(v3./f , 20*log10(abs(Fs)));
Ti = 0.25/v3;
z_1 = exp( - 1j*2*pi* f/v3 * 0.25);
% Cz3 = ( exp(Omega_1*Ti/2) - exp(- Omega_1*Ti/2)*z_1 )/(Omega_1*Ti);
% Gz3 = ( (1-z_1).^2 + (1 - z_1.^2)*a*Ti + (a*Ti + 2*a*Ti*z_1 + a*Ti*z_1.^2) ) ./ (Ti^2*a2b2);
biir = [1+a*Ti+aa*Ti*Ti ,-2*(1-aa*Ti*Ti), 1-a*Ti+aa*Ti*Ti];
Gz3 = (biir(1) + biir(2).*z_1 + biir(3).*z_1.^2) ./ ((aa+bb)*Ti*Ti);
hold on;semilogx(v3./f , 20*log10(abs(Gz3)));
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('\lambda /m');
ylabel('Mag /dB');
legend 50km/h 50km/h 120km/h 120km/h 350km/h 350km/h;
% legend 50km/h 120km/h 350km/h
grid on;

%% 单纯的画一遍去移变滤波器
figure1 = figure('Color',[1 1 1]);
semilogx( v1./f , 20*log10(abs(Gz1)));hold on;
semilogx( v2./f , 20*log10(abs(Gz2)));
semilogx( v3./f , 20*log10(abs(Gz3)));
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('\lambda /m')
ylabel('Mag /dB')
title('抗混叠滤波器+去移变滤波器')
legend 50km/h 120km/h 200km/h

%% 抗混叠+去移变
figure1 = figure('Color',[1 1 1]);
semilogx(v1./f , 20*log10(abs(Fs.*Gz1)));
hold on;
semilogx(v2./f , 20*log10(abs(Fs.*Gz2)));
semilogx(v3./f , 20*log10(abs(Fs.*Gz3)));
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('\lambda /m')
ylabel('Mag /dB')
title('抗混叠滤波器+去移变滤波器')
legend 50km/h 120km/h 200km/h

%% 新的研究： 车体振动
%%需要注意的是，这里就只到2Hz的位置就应该不画了
%%因为空间域的采样频率就是4Hz，假设车体速度为100m/s，那么采样率为100*4=400Hz
%%psi_max = f./v = 400Hz./v = 4Hz，所以还是4Hz
figure1 = figure('Color',[1 1 1]);
semilogx(f./v1 , 20*log10(abs(Fs.*Gz1)));
hold on;
semilogx(f./v2 , 20*log10(abs(Fs.*Gz2)));
semilogx(f./v3 , 20*log10(abs(Fs.*Gz3)));
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('\psi /Hz')
ylabel('Mag /dB')
title('抗混叠滤波器+去移变滤波器')
legend 50km/h 120km/h 200km/h
%%实际这里选点相当于低通滤波器，采样信号丢失算不算混叠了？因为高频的信号在低频采样时
%%(这里的采样是指空间域的选点)，因此认为还是存在信号的混叠。这里低通滤波器相当于可以
%%完成抗混叠的功能。

%%当然也可能并没有抗混叠的功能。

%%
%%这个抗混叠的功能还是存在的吧？

%% 验证一下这个抗混叠的功能是不是成立？
% figure;semilogx(f,20*log10(abs(Fs)));
%%为什么截止频率在1Hz左右，这个用法比较神奇

%%
% figure1 = figure('Color',[1 1 1]);
% semilogx(f, 20*log10(abs(Fs.*Gz1)));
% hold on;
% semilogx(f , 20*log10(abs(Fs.*Gz2)));
% semilogx(f , 20*log10(abs(Fs.*Gz3)));
% grid on;
% set(gca,'Fontname','Times New Roman','fontsize',14);
% xlabel('\psi /Hz')
% ylabel('Mag /dB')
% title('抗混叠滤波器+去移变滤波器')
% legend 50km/h 120km/h 200km/h


