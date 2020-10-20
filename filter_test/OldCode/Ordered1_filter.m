
% =========================================================================
%
%                  一阶加速度低通滤波器(抗混叠滤波器)与互补滤波器的问题探讨
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月29日
%   作者：
%--------------------------------------------------------------------------
%  功能： 1.首先对于传感器数据的处理，在抗混叠滤波器中的传感器速率是500
%        2.互补滤波器就不是了吧？
%       3. 
% 10.16
%       1. 信号的采样频率是500Hz，所以还是基于采样的
%       2. 滤波器参数随着速度变化，实际上就是随着TBS变化
%       3. 增加另一种表示
% 
% 
%--------------------------------------------------------------------------



%% 一阶数字抗混叠滤波器+补偿滤波器
clear all;
clc;
%----------一阶数字抗混叠滤波器------------
W1 = (10^5)/(2^17);
v1= 150/36;v2 = 200/9;v3 = 3500/36;
T1 = 0.25/v1;T2 = 0.25/v2;T3 = 0.25/v3;
%semilogx(f1,20*log10(abs(h1)));xlabel('频率（Hz）');ylabel('幅值(dB)');
figure();suptitle('一阶数字抗混叠+补偿');
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
xlabel('波长（m）');ylabel('幅值(dB)');


%%
%-----------一阶数字补偿滤波器----------------
T1 = 0.25/v1;
bc = [ 1+(W1*T1/2),  (-1)+(W1*T1/2)];
ac = [W1*T1 0];
[h1c f1c] = freqz(bc,ac,800000,500);
semilogx(v1./f1c,20*log10(abs(h1c)));hold on;

T2 = 0.25/v2;
bc = [ 1+(W1*T2/2),  (-1)+(W1*T2/2)];
ac = [W1*T2 0];
[h1c f1c] = freqz(bc,ac,800000,500);
semilogx(v2./f1c,20*log10(abs(h1c)),'g');hold on;

bc = [ 1+(W1*T3/2),  (-1)+(W1*T3/2)];
ac = [W1*T3 0];
[h1c f1c] = freqz(bc,ac,800000,500);
semilogx(v3./f1c,20*log10(abs(h1c)),'r');hold on;


%%
%-----------------一阶数字抗混叠+补偿滤波器-------
figure;
B1 = [2*W1*T1+(W1^2)*T1*T1 ,  (W1^2)*T1*T1-2*W1*T1];
A1 = [2*W1*T1+2*(W1^2)*T1*T1  ,  (-2)*W1*T1];
[H1 F1] = freqz(B1,A1,800000,500);
semilogx(v1./F1,20*log10(abs(H1)),'b');hold on;

B1 = [2*W1*T2+(W1^2)*T2*T2  , (W1^2)*T2*T2-2*W1*T2];
A1 = [2*W1*T2+2*(W1^2)*T2*T2  , (-2)*W1*T2];
[H1 F1] = freqz(B1,A1,800000,500);
semilogx(v2./F1,20*log10(abs(H1)),'g');hold on;

B1 = [2*W1*T3+(W1^2)*T3*T3 , (W1^2)*T3*T3-2*W1*T3];
A1 = [2*W1*T3+2*(W1^2)*T3*T3 , (-2)*W1*T3];
[H1 F1] = freqz(B1,A1,800000,500);
semilogx(v3./F1,20*log10(abs(H1)),'r');hold on;

xlabel('波长（m）');ylabel('幅值(dB)');


%% 这个不对也是很奇怪
figure; Cz = h1.*h1c;
semilogx(v3./F1,20*log10(abs(Cz)),'r');hold on;


%% 换一种表示
% f = 0.1:0.0001:100;
% fs = 500;
% z_1 = exp( - 1j*2*pi*f/fs);
% % s = (1-z^-1)/T，为什么这里的T用TBS取代？
% % 参数配置
% Omega_1 = (10^5)/(2^17);
% v1 = 150/3.6;
% T1 = 0.25/v1;
% b = [ Omega_1 * T1 0];
% a = [ 1 + Omega_1*T1 ,  (-1)];
% 
% % 画图
% Bz = Omega_1*T1./ ( 1 + Omega_1*T1 - z_1 );
% figure;semilogx(v1./f , 20*log10(Bz));
% 
% Cz = ( exp(Omega_1*T1/2) - exp(- Omega_1*T1/2)*z_1 )/(Omega_1*T1);
% hold on;semilogx(v1./f , 20*log10(Cz));




