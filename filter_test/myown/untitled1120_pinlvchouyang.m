

% =========================================================================
%
%                  测试
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月16日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.频率抽样法进行处理
%        2.
%       3. 
%--------------------------------------------------------------------------

%% 和下面的滤波器是一样的
%%修改自8.10这个例子的方法
%%这个滤波器的系数应该也可以变化的
clear all;
close all;
N = 40;
k = 0:6;
tp1 = exp(-1j*39/40*pi*k).*(7:-1:1);
tp2 = 0.388*exp(-1j*39/40*pi*7);
tp3 = zeros(1,length(8:32));
k = 34:39;
tp4 = 0.388*exp(-1j*39/40*pi*7);
tp5 = exp(-1j*39/40*pi*(N-k)).*(1:6);
Hk = [tp1,tp2,tp3,tp4,tp5];
hn = ifft(Hk);
figure1 = figure('Color',[1 1 1]);freqz(hn,1);
set(gca,'Fontname','Times New Roman','fontsize',14);

%% 验证滤波器的准确性
fs = 500;
htmp = fftshift(fft(hn,fs));
N = length(htmp);
x = (-N/2:N/2-1)/N*2;
figure;plot(x,20*log10(abs(htmp)));
list = htmp(find(x==0):find(x==0.3));
figure;plot(x(find(x==0):find(x==0.3)),abs(list));
%%确实与设想的差不多，0.1到0.3确实是线性递减的



%% 例子8.10
% clear all;
f = 1:0.01:250;
fs = 500;
omega = 2*pi*f/fs;
tmp = exp(- 1j.*omega);
% H = tmp.^(39/2)*(  );
t1 = sin(20*omega)/40./sin(omega/2);
t2 = zeros(1,length(omega));
for k =1:6
    t2 = t2 + sin(40*(omega/2-pi/40*k))./sin(omega/2-pi/40*k) + ...
        sin(40*(omega/2+pi/40*k))./sin(omega/2+pi/40*k);
end
t2 = t2/40;
t3 = 0.3/40*( sin(40*(omega/2-7*pi/40))./sin(omega/2-7*pi/40) + ...
    sin(40*(omega/2+7*pi/40))./sin(omega/2+7*pi/40));

H = tmp.^(39/2).*(t1+t2+t3);
figure1 = figure('Color',[1 1 1]);plot(f/fs*2,20*log10(abs(H)));grid on;
xlabel(' \omega (x\pi rad/sample)');
ylabel('Mag dB')
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);

%% 测试0.388，看上去都可以吧
% figure1 = figure('Color',[1 1 1]);
% for ii = 1:10
%     TEMP = 0.1*ii
%     %%
%     f = 1:0.01:250;
%     fs = 500;
%     omega = 2*pi*f/fs;
%     tmp = exp(- 1j.*omega);
%     % H = tmp.^(39/2)*(  );
%     t1 = sin(20*omega)/40./sin(omega/2);
%     t2 = zeros(1,length(omega));
%     for k =1:6
%         t2 = t2 + sin(40*(omega/2-pi/40*k))./sin(omega/2-pi/40*k) + ...
%             sin(40*(omega/2+pi/40*k))./sin(omega/2+pi/40*k);
%     end
%     t2 = t2/40;
%     t3 = TEMP/40*( sin(40*(omega/2-7*pi/40))./sin(omega/2-7*pi/40) + ...
%         sin(40*(omega/2+7*pi/40))./sin(omega/2+7*pi/40));
%     
%     H = tmp.^(39/2).*(t1+t2+t3);
%     if ii == 3
%         plot(f/fs*2,20*log10(abs(H)),'LineWidth',1);grid on;
% 
%     else 
%         plot(f/fs*2,20*log10(abs(H)));grid on;
%     end
%     hold on;
% end
% xlabel(' \omega (x\pi rad/sample)');
% ylabel('Mag dB')
% grid on;
% set(gca,'Fontname','Times New Roman','fontsize',14);
% legend 0.1 0.2 0.3
