
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
%  ���ܣ� 1.̽��һ�»��ֵ�����
%        2.
%        3. 
%--------------------------------------------------------------------------


%% ��������
clear all;
close all;
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
%%����Ĵ�������е�����
z_dot =  wavediff(1)/4;
z = 0;
for i = 1:length(z_acc)-1
    z_dot = z_dot + (z_acc(i) + z_acc(i+1))/2*dt*dt;
    z = z + z_dot;
    z_dot_save(i) = z_dot;
    z_save(i) = z;
end

% figure;plot(x,z_save);hold on;plot(x , longwave);legend 1 2
l2 = z_save - longwave;

%%
figure1 = figure('Color',[1 1 1]);plot(l1*1e3);hold on;plot(l2*1e3);
legend һ�׻��� ���λ���
set(gca,'Fontname','Times New Roman','fontsize',14);
%%�������ֵķ���ȷʵ��Խ������Ӱ�죬��ʱ��ȥ����


%% function
%%��������ȽϺ�ʱ
function out = cauDis(x1,x2)
    out = norm(x1 - x2);
end
function [out,out2 ] = cau(x_pos, f)
global L;

Dis = L;
Fx = [x_pos+Dis , f(x_pos+Dis)];
Bx = [x_pos , f(x_pos)];
out2 = cauDis(Fx,Bx);
out = (f(x_pos + Dis) - f(x_pos)) / Dis;
end
