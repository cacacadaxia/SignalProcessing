% =========================================================================
%
%                  �Ա������������ݣ����������
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 10��31��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.�����õ��Ǿ�������ݣ�������ʱû�в����ļ�����ʱ������ݵĶ��빤�����ҵ�
%               ��Ӧ�����ݣ��Թ�����ĳ����Աȷ���
%        2.С�����ݶ࣬��������١����������С���Ĳ��������0.25m��
% 11.24
%       3. 
%--------------------------------------------------------------------------

%% step1: ��ȡС������
clear;
close all;
data = load('gaodi_xiaoche_jj.txt');
gpproLTrolleyOri = data(:,2);
xSmallCar = data(:,1);
ll = xSmallCar(2:end)-xSmallCar(1:end-1);
b = load('filter_120m.mat');
len = (length(b.Num)-1)/2;
out = conv(b.Num,gpproLTrolleyOri);
gpproLTrolley = out(len+1:end-len);

%% step2: �ض�����
N = 5000;
gpproLTrolley = gpproLTrolley(81572 - N : 81572 + 1538 + N );
xSmallCar = xSmallCar(81572 - N : 81572 + 1538 + N);
gpproLTrolleyOriCut = gpproLTrolleyOri(81572 - N : 81572 + 1538 + N);

%% step3: �����Ƴ���0.65mת����Ϊ0.25m������ͬ��
gpproLTrolley25m = [];
for i = 1:floor(length(gpproLTrolley)*0.65*4)-4
    p1 = floor(i*0.25/0.65) + 1;
    p2 = mod(i*0.25 , 0.65);
    gpproLTrolley25m(i) = gpproLTrolley(p1) + p2/0.65*( gpproLTrolley(p1+1) - gpproLTrolley(p1) );
end
% figure;plot(0:0.25:(length(gpproLTrolley25m)-1)*0.25 , gpproLTrolley25m,'-r');
% hold on;plot(0:0.65:(length(gpproLTrolley)-1)*0.65 , gpproLTrolley,'--k');

%% step4: ��ȡ��쳵����
 data2 = load('gaodi_guijian_jj.txt');
 gpproGjOri = data2(:,6);
 xGj = data2(:,1);
% figure1 = figure('Color',[1 1 1]);
% plot(xGj,gpproGjOri);
% %%������ݶ����ˣ�����excel�е������Ϣ
% hold on;plot(xSmallCar , gpproLTrolley);
% set(gca,'Fontname','Times New Roman','fontsize',14);

%% step5: ����������ҵ������
out = conv( gpproLTrolley25m , gpproGjOri(end:-1:1));
[~,index] = max(out);
index = 17875;
index = index - length(gpproGjOri);
figure1 = figure('Color',[1 1 1]);
plot( gpproLTrolley25m);hold on;plot([zeros(1,index) , gpproGjOri'])
xlabel('������ 0.25m');ylabel('�ߵ� mm');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
legend С�� ���


%% step6: ��������
gpproLTrolleyOriCut = gpproLTrolleyOri;
gpproLTrolley25m = [];
for i = 1:floor(length(gpproLTrolleyOriCut)*0.65*4)-4
    p1 = floor(i*0.25/0.65) + 1;
    p2 = mod(i*0.25 , 0.65);
    gpproLTrolley25m(i) = gpproLTrolleyOriCut(p1) + p2/0.65*( gpproLTrolleyOriCut(p1+1) - gpproLTrolleyOriCut(p1) );
end
b = load('filter_250m.mat');
len = (length(b.Num)-1)/2;
out = conv(b.Num,gpproLTrolley25m);
gpproLTrolley_1000m = out(len+1:end-len);
% figure1 = figure('Color',[1 1 1]);plot(gpproLTrolleyOriCut+196);hold on ;plot(gpproLTrolley);
% legend ԭʼ���� �˲���;
% set(gca,'Fontname','Times New Roman','fontsize',14);


%% step7: ������ͬ�����Ĳ��ζԱ�
b = load('filter_120m.mat');
len = (length(b.Num)-1)/2;
out = conv(b.Num,gpproLTrolley25m);
gpproLTrolley_120m = out(len+1:end-len);
figure1 = figure('Color',[1 1 1]);plot(gpproLTrolley_120m);hold on ;plot(gpproLTrolley_1000m);
legend 120m 250m;
set(gca,'Fontname','Times New Roman','fontsize',14);


%% step8������Ƶ��
% figure;plot(gpproLTrolleyOri);
% plot_mag(gpproLTrolleyOri,'���С������');
% plot_mag(gpproLTrolley_1000m,'���С������');
% plot_mag(gpproLTrolley_120m,'���С������','hold');
% legend 250m 120m
%% ��ͷ���ٶ�
% wave1 = gpproLTrolleyOri./1e3;
% N = 5*4;%%5m
% pitch_t_1 = atan( (wave1(N+1:end) - wave1(1:end-N))./0.25./N )/pi*180;
% N = 10*4;%%10m
% pitch_t_2 = atan( (wave1(N+1:end) - wave1(1:end-N))./0.25./N )/pi*180;
% figure1 = figure('Color',[1 1 1]);
% plot(pitch_t_1);hold on;plot(pitch_t_2);legend L=5m L=10m
% xlabel('������ (0.25m)');ylabel('�Ƕ� (deg)');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% 
% %%�ɼ�����ĵ�ͷ�Ƕȶ�ͦС�ģ��������ڲ����Ƚϴ�������

%%
% plot_mag(pitch_t_1,'�Ƕ�Ƶ��');
% plot_mag2(pitch_t_1,'�Ƕ�Ƶ��');


%%
function plot_mag(signal_data , tit , varargin)
if (nargin == 3)
    mode = varargin{1};
    if mode == 'hold'
        hold on;
    end
else
    figure1 = figure('Color',[1 1 1]);
end

%% ������һ��
% E = sum(abs(signal_data).^2);
% signal_data = signal_data./sqrt(E);

%%
fs = 4;     %% 0.25mΪһ���������
signal_data = signal_data./1e3;
N = length(signal_data);
x = (1:N/2+1)/N*fs;
x = 1./x;
tp = abs(fftshift(fft(signal_data)));
tp = 20*log10(tp);
semilogx( x , tp((length(tp)/2):end)   );
xlabel('\lambda m')
ylabel('Mag dB')
set(gca,'Fontname','Times New Roman','fontsize',14);
title(tit);
grid on;
end

function plot_mag2(signal_data , tit , varargin)
if (nargin == 3)
    mode = varargin{1};
    if mode == 'hold'
        hold on;
    end
else
    figure1 = figure('Color',[1 1 1]);
end

%% ������һ��
% E = sum(abs(signal_data).^2);
% signal_data = signal_data./sqrt(E);

%%
fs = 4;     %% 0.25mΪһ���������
signal_data = signal_data./1e3;
N = length(signal_data);
x = (1:N/2+1)/N*fs;
x = 1./x;
tp = abs(fftshift(fft(signal_data)));
tp = 20*log10(tp);
plot( x , tp((length(tp)/2):end)   );
xlabel('\lambda m')
ylabel('Mag dB')
set(gca,'Fontname','Times New Roman','fontsize',14);
title(tit);
grid on;
end


