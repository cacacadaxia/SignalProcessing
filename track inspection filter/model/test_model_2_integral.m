
% =========================================================================
%
%                  ���������
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 10��20��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.̽��һ�»��ֵ����⣬�Ƚϸ��ֻ��ֵĽ����
% ʵ���ϸ�����ɢ�Ļ��ֶ�����д�ɴ��ݺ�������ʽȥ��⣬���յ����
%        2.
%        3. 
%--------------------------------------------------------------------------

%% ��������
clear all;
close all;
% clear ydot ydot_1 y_1;
global L;
L = 3;
del_x = 0.25;               %%��λΪm
waveLen = 300;              %%��λΪm
x = 0:del_x:waveLen*10;     %%�ռ�������
waveMag = 40;               %%mm �����ķ�ֵ
longwave = waveMag*1e-3*sin(2*pi/waveLen*x);
% figure;plot(x,longwave);xlabel('������� /0.25m');ylabel('��ֵ /mm');

%% Ϊʲôȥ�ҳ����������أ�
%%��Ϊ�������ߵ����߾������Ź��ܼ�������ĵ�ͷ���ٶ�
wavediff = waveMag*1e-3*2*pi/waveLen*cos(2*pi/waveLen* x );

%% ʱ������
v = 100/3.6;                                             %%m/s
dt = del_x/v;                                       %%����ʱ����
t = 0:dt:(length(x)-1)*dt;
T = dt*waveLen/del_x;                               %%ʱ���������
z_true = waveMag*1e-3*sin(2*pi/T*t);                %%ʱ������
z_v = waveMag*1e-3*2*pi/T*cos(2*pi/T*t);            %%�ٶ�����
z_acc = -waveMag*1e-3*2*pi/T*2*pi/T*sin(2*pi/T*t);  %%���ٶȼƵ�ֵ


%% ���ٶȼ�
z_dot =  wavediff(1)/4;
z = 0;
for i = 1:length(z_acc)
    z_dot = z_dot + z_acc(i)*dt*dt;
    z = z + z_dot;
    z_dot_save(i) = z_dot;
    z_save(i) = z;
end
% figure;plot(x,z_save);hold on;plot(x , longwave);legend 1 2
l1 = z_save-longwave;

%% ���λ���
% %%����Ĵ�������е�����
% z_dot =  wavediff(1)/4;
% z = 0;
% for i = 1:length(z_acc)-1
%     z_dot = z_dot + (z_acc(i) + z_acc(i+1))/2*dt*dt;
%     z = z + z_dot;
%     z_dot_save(i) = z_dot;
%     z_save(i) = z;
% end
% 
% % figure;plot(x,z_save);hold on;plot(x , longwave);legend 1 2
% l2 = z_save - longwave;

%% ������
sigma_acc = 1e-5 * 9.8  / sqrt(dt);
z_acc = z_acc + randn(1,length(t)) * sigma_acc;


%%
% x_dot_0 =  2.094395102393196e-04;
x_dot_0 = wavediff(1)/4;
for i = 1:length(z_acc)
    z_save(i) = double_integral(z_acc(i)*dt*dt , x_dot_0);
end
l3 = z_save - longwave;



% figure;plot(z_save);hold on;plot(longwave(2:end));legend 1 2
l4 = z_save(1:end-1) - longwave(2:end);
figure1 = figure('Color',[1 1 1]);
plot(l3*1e3);hold on;plot(l4*1e3);legend ����ǰ ������;
set(gca,'Fontname','Times New Roman','fontsize',14);
ylabel('��� /mm');xlabel('��� /0.25m');

%%
% figure1 = figure('Color',[1 1 1]);plot(l1*1e3);hold on;plot(l2*1e3);plot(l3*1e3,'LineWidth',2);
% legend һ�׻��� ���λ���
% set(gca,'Fontname','Times New Roman','fontsize',14);
%%�������ֵķ���ȷʵ��Խ������Ӱ�죬��ʱ��ȥ����

%% ���λ���
function y_k = double_integral(x_k , x_dot_0)
persistent ydot ydot_1 y_1;
if isempty(ydot)
    ydot = 0;
    ydot_1 = x_dot_0;%%��ֵ����Ҫ
    y_1 = 0;
end
ydot = ydot_1 + x_k;
y_k = y_1 + ydot;
%%
y_1 = y_k;
ydot_1 = ydot;
end
    
