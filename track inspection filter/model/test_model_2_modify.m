
% =========================================================================
%
%                  ���������
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 11��4��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.���ϸ��ַ�����������ʱ������
%        2.
%        3. 
%--------------------------------------------------------------------------


%% ��������
clear all;
close all;
global L;
L = 3;
del_x = 0.25;               %%��λΪm
waveLen = 40;              %%��λΪm
waveLen = 10;              %%��λΪm
% waveLen = 20;
x = 0:del_x:waveLen*10;     %%�ռ�������
waveMag = 40;               %%mm �����ķ�ֵ
waveMag = 10;               %%mm �����ķ�ֵ
longwave = waveMag*1e-3*sin(2*pi/waveLen*x);
% figure;plot(x,longwave);xlabel('������� /0.25m');ylabel('��ֵ /mm');

%% Ϊʲôȥ�ҳ����������أ�
%%��Ϊ�������ߵ����߾������Ź��ܼ�������ĵ�ͷ���ٶ�
% wavediff = waveMag*1e-3*2*pi/waveLen*cos(2*pi/waveLen* x );
wavediff = (longwave(2:end) - longwave(1:end-1))./del_x;

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
l1 = z_save(1:end-1) - longwave(2:end);
l1 = [l1 , 0];

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
l2 = z_save(1:end) - longwave(1:end) - 1.046956521166820e-04;

figure1 = figure('Color',[1 1 1]);plot(x,l2*1e3,'r');hold on;plot(x,l1*1e3,'k--');legend ���ݼ�����ײ�� ���ٶȼƼ�����ײ��;
xlabel('������ /0.25m');ylabel('��� /mm')
set(gca,'Fontname','Times New Roman','fontsize',16);
%%����������˵������Ӧ������ȵģ�û��̫��Ĳ���
%%������������Ĳ����ǻ��ִ�������
%%���������Ϊ�����궨�Ĳο�����

figure;plot(z_save(1:end-1)  );hold on;plot(longwave(2:end));

%% �������������
f = @(x)(waveMag*1e-3*sin(2*pi/waveLen*x));%%����
for i = 1:length(x)
    [pitch2(i) , out2(i)] = cau(x(i),f);
end
% figure1 = figure('Color',[1 1 1]);plot(pitch2);hold on;plot(pitch(7:end));legend 1 2;
%%��˵����ʲô����˵����"����3m�����Ĺ�����н���ʵ�ʵ��н����"�ǳ�С�����Ժ��Բ���
%%ֻ��������֮������ʱ��Ϊʲô����Ϊʵ����Ҳ�����ߣ���(z1+z2)/2������
%%Ϊʲô��ZM�������أ�

%% �����(z1+z2)/2�Ĳ�ƽ˳���
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
z_dot = 2.093017317643780e-04;%%����Ϊ300m
% z_dot = 2.093017317643780e-04*15.3;
% z_dot = 2.093017317643780e-04*9.5;
z_dot = 3.366695773512057e-04*2;%%20m
z_dot = 3.366695773512057e-04*1.1236;%%40m
z_dot = 3.366695773512057e-04*2.3536;%%10m

z = 0;
for i = 1:length(wy2)
    z_dot = z_dot + wy2(i) * dt * del_x;
    z = z + z_dot;
    z_dot_save(i) = z_dot;
    z_save(i) = z;%%�����z_save�Ͳ����Ǹ������%%����Ľ���ܵ�������ֵ�ı仯Ӱ��Ƚϴ�
    %%���������Ǻͼ��ٶȼƵķ����յ��ĸ��žͺ�С
    %%����Ҳ������⣬��������Ƕ��Ǻͷ�ֵ��صġ�
    %%���������ǣ��������40mm�ķ�ֵ���У�������һ�Ŷ�����
end
det2 = z_save * 2;
y1 = (det2 - det1)/2;
l3 = y1(1:end-1) - longwave(2:end) - 1.046956521166820e-04;
l3= [l3 , l3(end)];
% figure1 = figure('Color',[1 1 1]);plot(x,(y1 - longwave)*1e3+1.05);
% xlabel('������ /0.25m');ylabel('��� /mm');
% set(gca,'Fontname','Times New Roman','fontsize',16);
%%�����з�ֵ�Ĳ�������Ǵ�������ģ�������Ҫ���²���


%% z1+z2/2==? pitch
figure1 = figure('Color',[1 1 1]);
plot((z_save+0.005068982715618)*1e3);
hold on;
l0 = f_tp(x);
plot(l0*1e3)
xlabel('������ (0.25m)');ylabel('��λ m');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
title('����');
legend result (z1+z2)/2
cauMag(l0)/2*1e3/waveMag
cauMag(z_save)/2*1e3/waveMag


%% �Ա����еķ���
% y1 = (y1(2:end)+y1(1:end-1))/2;
% y1 = [0,y1];
% figure1 = figure('Color',[1 1 1]);plot( x , l3*1e3 + 1.25);hold on;
% plot(x,l2*1e3,'r');hold on;plot(x,l1*1e3,'k--');legend ���������  ���ݼ�����ײ�� ���ٶȼƼ�����ײ��;
% xlabel('������ /0.25m');ylabel('��� /mm');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;

%%��������֮�󣬷����������һ�ײ����͵Ľ����һ����
%%�����Ա���Ȼ������һ����������������ܸо��е����⣬��ô���ܺͳ�����ϵ���أ�

%% �µķ����Ľ����ʾ
% figure1 = figure('Color',[1 1 1]);
% plot(y1*1e3+1.25);hold on;plot(longwave*1e3)
% xlabel('������ /0.25m');ylabel('�ߵͲ�ƽ˳ /mm')
% set(gca,'Fontname','Times New Roman','fontsize',14);
% grid on;
% legend result ref;
% title('�·����Ľ��');

function out = cauMag(in)
max_ = 155;
min_ = 115;
out = max(in(min_:max_))-min(in(min_:max_));
end



%% function
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
