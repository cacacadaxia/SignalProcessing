
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
%  ���ܣ� 1.
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

%% һ�ײ����͵Ĺ̶�������£�
% z = 0;
% for i = 1:length(z_v)
%     z = z + z_v(i)*dt;
%     z_save(i) = z;
% end
% % figure;plot(x,z_save);hold on;plot(x , longwave);legend 1 2
% l1 = z_save-longwave;
%%����˵������һ�ײ��ȥ��ͻ��й̶����������1/50����


%% ������
% wavediff_k_1 %%w*del_t
pitch = atan(wavediff);
wy = ( pitch(2:end) - pitch(1:end-1) )/dt;  %%�����ǲ����������
wy = [0,wy ];                               %%���Ӱ��ܶ�

z_dot = wavediff(1)/4;
z = 0;
for i = 1:length(wy)
    z_dot = z_dot + wy(i) * dt * del_x;
    z = z + z_dot;
    z_dot_save(i) = z_dot;
    z_save(i) = z;
end
l2 = z_save - longwave;
% figure1 = figure('Color',[1 1 1]);plot(x,l2*1e3,'r');hold on;plot(x,l1*1e3,'k--');legend ���ݼ�����ײ�� ���ٶȼƼ�����ײ��;
xlabel('������ /0.25m');ylabel('��� /mm')
set(gca,'Fontname','Times New Roman','fontsize',16);
%%����������˵������Ӧ������ȵģ�û��̫��Ĳ���
%%������������Ĳ����ǻ��ִ�������
%%���������Ϊ�����궨�Ĳο�����

%% �������������
f = @(x)(waveMag*1e-3*sin(2*pi/waveLen*x));

for i = 1:length(x)
    [pitch2(i) , out2(i)] = cau(x(i),f);
end
% figure1 = figure('Color',[1 1 1]);plot(pitch2);hold on;plot(pitch(7:end));legend 1 2;
%%��˵����ʲô����˵����"����3m�����Ĺ�����н���ʵ�ʵ��н����"�ǳ�С�����Ժ��Բ���
%%ֻ��������֮������ʱ��Ϊʲô����Ϊʵ����Ҳ�����ߣ���(z1+z2)/2������




%% 
% (z1+z2)/2
f_tp = @(x)(( waveMag*1e-3*sin(2*pi/waveLen*x) + waveMag*1e-3*sin(2*pi/waveLen*(x+L)) )/2);
wavediff_tp = (waveMag*1e-3*2*pi/waveLen*cos(2*pi/waveLen*x) + waveMag*1e-3*2*pi/waveLen*cos(2*pi/waveLen*(x+3)))/2;
% figure;plot((atan(wavediff_tp) - pitch2)/pi*180);
%%����˵����ʵ���ϲ����(z1+z2)/2�Ľ��

%% �����µķ���
det1 = sin(pitch2)*L;

num = 0;
% det1 = [zeros(1,num),det1(1:end-num)];
% det1 = [det1(1+num:end),zeros(1,num)];
% wavediff_k_1 %%w*del_t

wy2 = ( pitch2(2:end) - pitch2(1:end-1) )/dt;%%�����ǲ����������
wy2 = [0 , wy2];

% z_dot = wavediff_tp(1)/4;%%ֻ�����������
z_dot = 2.093017317643780e-04;
z = 0;
for i = 1:length(wy2)
    z_dot = z_dot + wy2(i) * dt * del_x;
    z = z + z_dot;
    z_dot_save(i) = z_dot;
    z_save(i) = z;
end
det2 = z_save * 2;
y1 = (det2 - det1)/2;
l3 = y1 - longwave;
% figure1 = figure('Color',[1 1 1]);plot(x,(y1 - longwave)*1e3+1.05);
xlabel('������ /0.25m');ylabel('��� /mm');
set(gca,'Fontname','Times New Roman','fontsize',16);
%%�����з�ֵ�Ĳ�������Ǵ�������ģ�������Ҫ���²���


%% �Ա����еķ���
% y1 = (y1(2:end)+y1(1:end-1))/2;
% y1 = [0,y1];
figure1 = figure('Color',[1 1 1]);plot( x , l3*1e3 + 1.25);
hold on;
plot(x,l2*1e3,'r');hold on;plot(x,l1*1e3,'k--');legend ���������  ���ݼ�����ײ�� ���ٶȼƼ�����ײ��;
xlabel('������ /0.25m');ylabel('��� /mm')
set(gca,'Fontname','Times New Roman','fontsize',16);

%%��������֮�󣬷����������һ�ײ����͵Ľ����һ����
%%�����Ա���Ȼ������һ����������������ܸо��е����⣬��ô���ܺͳ�����ϵ���أ�




%%
figure1 = figure('Color',[1 1 1]);plot(x(1:end-6),(z_save(1:end-6) - longwave(7:end))*1e3);
figure;plot(z_save);hold on;plot(longwave);legend 1 2


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
