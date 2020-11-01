
% =========================================================================
%
%                  ���ֹ�������㷨����
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��23��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.�ο�������˲��� / ����
%        2.
%         %% ��ƫ�˲���
%         %% ��ƫ�˲������ܵ��ٶȵ�Ӱ��
%         %% Ϊʲô����Ϊ����˲������������ٶȱ仯��ֹƵ�ʵģ������Ǻ�ǰ���˲���һ�µ�
%         %% ������ͼ�����յ�Ӱ��
% 
% 
% 2.1: ��Ϊǰ�˲ɼ��źŵ������ԣ�Ҳ�ǰ���ʱ�����вɼ�Ȼ����ѡ��
%        3. sin��exp�����ֱ�ʾ����
% 
% 
%--------------------------------------------------------------------------


%% ��϶����߶���һ���ģ��ѵ������ʲ�Ӧ����֮�仯��
%%��������߼�
close all;
clear all;
wd = 0.001;
b = [1-wd,-1+wd];
a = [1,-1+wd];
% b = [1,2,3];
% a = [4,5];

v = 100/3.6;
t1 = 0.25/v;
[h,f] = freqz(b,a,10000,1/t1);
lamda = v./f;
figure;
semilogx(lamda,20*log10(abs(h)));

v = 200/3.6;
t1 = 0.25/v;
[h,f] = freqz(b,a,10000,1/t1);
lamda = v./f;
hold on;
semilogx(lamda,20*log10(abs(h)));
legend 1 2

v = 300/3.6;
t1 = 0.25/v;
[h,f] = freqz(b,a,10000,1/t1);
lamda = v./f;
hold on;
semilogx(lamda,20*log10(abs(h)));
% legend 1 2 3

%% �����ʲ���һ����������������ͬһ������
lamda = 1:0.001:10000;
psi = 1./lamda;
wd = 0.001;
Ls = 0.25;
Rz = (1-wd)*(1- exp(-j*2*pi*psi*0.25))./( 1-(1-wd)*exp(-j*2*pi*psi*0.25) );
hold on;semilogx(lamda,20*log10(Rz),'LineWidth',1);

%% ����ɶ��˼��û�����ף��������ҵ������ʽ��
tmp = 2*(1-wd)^2*(1-cos(2*pi.*psi*0.25))./( 1+ (1-wd)^2 -2*(1-wd)*cos(2*pi.*psi*0.25) );
Rz2 = 10*log10(tmp);
%%bode
figure;semilogx(1./psi , 20*log10(Rz),'r--','LineWidth',1);
hold on;semilogx(1./psi , Rz2,'LineWidth',1);

%% ģ����һ��freqz�Ĺ���
% %%��0��25�޹�
% psi = 0:0.001:900;
% fs = 1000;
% wd = 0.001;
% b = [1-wd,-1+wd];
% a = [1,-1+wd];
% z_1 = exp(-j*2*pi*psi/fs);
% Rz = (1-wd)*(1 - z_1)./( 1-(1-wd) * z_1 );
% figure;semilogx(psi,20*log10(Rz),'LineWidth',1);
% [h,f] = freqz(b,a,1e5,fs);
% hold on;semilogx(f,20*log10(h),'--r','LineWidth',2)


%% �����˲���(ͨ��fdatool��Ƶ��˲���)
b = load('filter3.mat');
[h,f] = freqz(b.Num,[zeros(1,length(b.Num)-1),1],10000,4);
figure1 = figure('Color',[1 1 1]);semilogx(1./f,20*log10(h));grid on;xlabel('���� /m');
ylabel('��ֵ /dB')
set(gca,'Fontname','Times New Roman','fontsize',16);

% figure;semilogx(1./f,angle(h));grid on;





