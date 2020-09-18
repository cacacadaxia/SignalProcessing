
% =========================================================================
%
%                  ���ֹ�������㷨����
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��16��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.���ȸ�����Ϊ�򵥵Ĺ��򲿷�
%        2.�����˲����Ĳ��֣������жԱȡ���Ҫ�ǹ���
%        3. 
%--------------------------------------------------------------------------

% load_txt;
close all;
clear all;

load_txt;
size(wave_out);
N = length(fmctrl_data);
x = 0:0.25:0.25*(N-1);
x = x/1000;

%% ���ĶԱ�
chaogao = wave_out(:,6);

% ****************�����趨*********************
delay = 418;%%һֱ���������
%%������������ɶ��˼������������׼ȷ��
G_par = 3.8259e-14*141500.0;
ht = 3.90398e-05*141500.0*0.268;

tbs = fmctrl_data(:,end);
tbs_s = tbs/1e5;
% ***************step1 ģ�ʹ*************************
rou_l = fmctrl_data(:,9);
rou_r = fmctrl_data(:,11);

% rou_l = rou_l/(129.01);
% rou_r = rou_r/(129.01);

for i = 3:length(rou_l)
    rou_l_dot2(i,1) = rou_l(i) - 2*rou_l(i-1) + rou_l(i-2);
    rou_r_dot2(i,1) = rou_r(i) - 2*rou_r(i-1) + rou_r(i-2);
end


% ******************step2 ���ٶȼƵ��˲� **********************************
ay = fmctrl_data(:,5);
% ay = ay/(32768/(2*9.83));

% ---------------------- �����˲��� --------------------
for i = 1:length(ay)
    ay_Gz(i,1) = G(ay(i) , tbs(i));
end

%% Ƶ�׹۲�
% plot_mag(ay,'�˲�ǰ')
% plot_mag(ay_Fz,'�˲���')
% plot_mag(wave_out(:,3),'�����')
% 
%% ����
% sita_b = sita_b/3276.8/180*pi;
sita_b = aln(:,4);
for i = 3:length(sita_b)
    sita_b_dot2(i,1) = sita_b(i) - 2*sita_b(i-1) + sita_b(i-2);
end

camo = ay_Gz - G_par .* sita_b .* tbs.^2 + ht*sita_b_dot2;
camo = -camo;

yL_dot2 = camo - rou_l_dot2;
yR_dot2 = camo + rou_r_dot2;


x_dot = 0;
x = 0;

for i = 3:length(yL_dot2)           %%�򵥻��֣��϶��ǲ��Ե�
    x_dot = x_dot + yL_dot2(i)*tbs(i);
    x = x + x_dot*tbs(i);
    yL(i,1) = x;
end


% Ŀǰֻ�л����Ǵ��������

%% ����Ա�

figure;plot(camo(30:end));hold on;
plot(aln(:,1));
legend 1 2

figure;plot(yL);



%% ����Ĳ���Ƶ��
aln_l = wave_out(:,3);
aln_r = wave_out(:,4);

function out = F(x,tbs)
%% �˲����趨
%% ����x[3]������ͬ��
% x(3)=x_n
% x(2)=x_n-1
% x(1)=x_n-2
%% ����˲���������ʱ����ʮ����
persistent y;
if isempty(y)
    y = zeros(3,1);
end
y(3) = ( y(2)*(2*2^28+2^14*tbs) - y(1)*2^28+tbs^2*x  )/(2^28 + 2^14*tbs + tbs^2);
%% ����temp��
y(1) = y(2);
y(2) = y(3);
out = y(3);
end


function out = R(x_k,tbs)
wd = 0.001;

persistent y x;
if isempty(y)
    y = zeros(2,1);
    x = zeros(2,1);
end
x(2) = x_k;
x_dot = x(2)-x(1);
y(2) = (1 - wd) * (x_dot + y(1));

%% ����temp��
x(1) = x(2);
y(1) = y(2);
out = y(2);

end

function out = G(x_k,tbs_k)
persistent x y1 tbs;
if isempty(x)
%     y = zeros(2,1);
    x = zeros(3,1);
    y1 = zeros(2,1);
    tbs = zeros(2,1);
end
%% ����kʱ��
x(3) = x_k;
tbs(2) = tbs_k;

%% ��ɢ�˲�
x_dot = x(3) - x(2);
x_dot_2 = x(3) - 2*x(2) + x(1);
y1(2) = (x(3) + x(2)) *tbs(2)/2^15 + x_dot;
y = (y1(2)+y1(1))* (tbs(2) + tbs(1))/2^16 + x_dot_2;


%% ����temp��
y1(1) = y1(2);
x(1) = x(2);
x(2) = x(3);
tbs(1) = tbs(2);
out = y;

end

function plot_mag(signal_data , tit)
figure;
fs = 4;     %% 0.25mΪһ���������
N = length(signal_data);
x = (1:N/2+1)/N*fs;
x = 1./x;
tp = abs(fftshift(fft(signal_data)));
tp = 20*log10(tp);
semilogx( x , tp((length(tp)/2):end)   );
xlabel('\lambda m')
ylabel('Mag dB')
set(gca,'Fontname','Times New Roman','fontsize',16);
title(tit);

end










