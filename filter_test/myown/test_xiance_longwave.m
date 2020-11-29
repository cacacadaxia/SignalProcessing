

% =========================================================================
%
%                  弦测法的讨论
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 10月23日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.弦测法设计
%        2.
%       3. 
%--------------------------------------------------------------------------

clear all;
close all;

%% 10m弦测
psi = 1/200:0.001:1/4;
lambda = 1./psi;
% lambda = 0.5:0.1:200;%%最小波长为0.5m
L = 10;
HlambdaMag = 1 - cos(pi./lambda*L);
Hfunc = @(lambda)(1 *(1 - cos(pi./lambda*L)));
% figure1 = figure('Color',[1 1 1]);loglog((psi),(HlambdaMag));xlabel('频率 (1/m) ');ylabel('M(x)/f(x)')
figure1 = figure('Color',[1 1 1]);plot((lambda),(HlambdaMag));xlabel('波长 m ');ylabel('幅值响应 dB');set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
figure1 = figure('Color',[1 1 1]);plot((psi),(HlambdaMag));xlabel('频率 (1/m) ');ylabel('M(x)/f(x)');grid on;
figure1 = figure('Color',[1 1 1]);plot((psi),(1./HlambdaMag));xlabel('频率 (1/m) ');ylabel('M(x)/f(x)');grid on;
Ifunc = @(lambda)(1 ./ (1 - cos(pi./lambda*L)));

%% 逆傅立叶变换
%%20--100
psi = 1/200:0.001:1/10;
lambda = 1./psi;
HlambdaMag = 1 - cos(pi./lambda*L);
Hz = 1./HlambdaMag;
figure1 = figure('Color',[1 1 1]);
plot(psi,(Hz));
xlabel('\psi 1/m');
ylabel('Mag ');
hold on;
plot(psi,1/100*psi.^(-1.65));%%拟合函数是这个
legend 需要设计的滤波器的幅值 拟合的幅值
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
%% 只是单纯的通带滤波器
% psi_tmp = 0.01:0.005:0.1;
% % filter_coeff1 = 1/100*psi_tmp.^(-1.65);
% filter_coeff = 1./(1 - cos(pi*psi_tmp*L));
% N = 800;
% tp1 = zeros(1,2);
% k = 2:20;
% tp2 = exp(-1j*(N-1)/N*pi*k);
% tp3 = 0.388*exp(-1j*(N-1)/N*pi*21);
% k = N-21;
% tp4 = zeros(1,length(22:k-1));
% tp5 = 0.388*exp(-1j*(N-1)/N*pi*k);
% k = N-20:N-2;
% tp6 = exp(-1j*(N-1)/N*pi*k);
% tp7 = [0];
% Hk = [tp1, tp2,tp3,tp4,tp5,tp6,tp7];
% hn = ifft(Hk);
% fs = 4;
% [H,f] = freqz(hn,1,1e4,fs);
% 
% figure1 = figure('Color',[1 1 1]);plot(f,20*log10(abs(H)));grid on;
% xlabel(' f Hz');
% ylabel('Mag dB');
% hold on;
% psi_tmp = 0.01:0.005:0.1;
% lambda = 1./psi_tmp;
% HlambdaMag = 1 - cos(pi./lambda*L);
% Hz = 1./HlambdaMag;
% plot(psi_tmp , 20*log10(Hz),'LineWidth',1);
% grid on;
% set(gca,'Fontname','Times New Roman','fontsize',14);


%% 基于频率抽样法的FIR设计
% psi_tmp = 0.01:0.005:0.1;
% % filter_coeff1 = 1/100*psi_tmp.^(-1.65);
% filter_coeff = 1./(1 - cos(pi*psi_tmp*L));
% N = 800;
% tp1 = zeros(1,2);
% k = 2:20;
% tp2 = exp(-1j*(N-1)/N*pi*k).*filter_coeff;
% tp3 = 0.388*exp(-1j*(N-1)/N*pi*21);
% k = N-21;
% tp4 = zeros(1,length(22:k-1));
% tp5 = 0.388*exp(-1j*(N-1)/N*pi*k);
% k = N-20:N-2;
% tp6 = exp(-1j*(N-1)/N*pi*k).*filter_coeff(end:-1:1);
% tp7 = [0];
% Hk = [tp1, tp2,tp3,tp4,tp5,tp6,tp7];
% hn = ifft(Hk);
% fs = 4;
% [H,f] = freqz(hn,1,1e4,fs);
% 
% figure1 = figure('Color',[1 1 1]);plot(f,20*log10(abs(H)));grid on;
% xlabel(' f Hz');
% ylabel('Mag dB');
% hold on;
% psi_tmp = 0.01:0.005:0.1;
% lambda = 1./psi_tmp;
% HlambdaMag = 1 - cos(pi./lambda*L);
% Hz = 1./HlambdaMag;
% plot(psi_tmp , 20*log10(Hz),'LineWidth',1);
% grid on;
% set(gca,'Fontname','Times New Roman','fontsize',14);
% 
% %%对比发现其结果没有比较大的问题
% % save('xiancefilter.mat','hn')


%% 论文中的方法
%%但是两者对于弦测方法的定义并不一致
%%论文中的是两点差分法，所以暂时先认为上面是对的。
del_zeta = 0.005;
psi_tmp = 0.005:0.005:0.1;
a = 2;
filter_coeff = 1./( 2*sin(2*pi*psi_tmp) );
N = 800;
tp1 = zeros(1,1);
k = 1:20;
tp2 = exp(-1j*(N-1)/N*pi*k).*filter_coeff;
tp3 = 0.388*exp(-1j*(N-1)/N*pi*21);
k = N - 21;
tp4 = zeros(1,length(22:k-1));
tp5 = 0.388*exp(-1j*(N-1)/N*pi*k);
k = N-20:N-1;
tp6 = filter_coeff(end:-1:1).* exp(-1j*(N-1)/N*pi*k);
tp7 = [0];
Hk = [tp1, tp2,tp3,tp4,tp5,tp6];
hn = ifft(Hk);
fs = 4;
[H,f] = freqz(hn,1,1e4,fs);

figure1 = figure('Color',[1 1 1]);plot(f,20*log10(abs(H)));grid on;
xlabel(' f Hz');
ylabel('Mag dB');
hold on;
psi_tmp = 0.01:0.005:0.1;
lambda = 1./psi_tmp;
Hz = 1./( 2*sin(2*pi*psi_tmp) );
plot(psi_tmp , 20*log10(Hz),'LineWidth',1);
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);

%%
figure1 = figure('Color',[1 1 1]);plot(1./f,20*log10(abs(H)));grid on;
xlabel(' 波长 m');
ylabel('Mag dB');
hold on;
psi_tmp = 0.005:0.005:0.1;
Hz = 1./( 2*sin(2*pi*psi_tmp) );
plot(1./psi_tmp , 20*log10(Hz),'LineWidth',1);
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);






