% =========================================================================
%
%                  ������ⷽ������
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 10��28��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.����һ�²�ƽ˳�������ø���һЩ
%        2.������Ҫ�յ���������С��Ӱ�죬���Ǽ��ٶȼ��յ��ٶ�Ӱ��Ƚϴ�(JC-21���ֻ����Ե���ϵͳ����ָ��JC-21���ֻ����Ե���ϵͳ����ָ��.docx)
%        3. 
%--------------------------------------------------------------------------

clear all;
close all;
% lamda1 = 8;
del_x = 0.25;
ftime = @(t,v)(1*(2*1e-3*sin(2*pi/8.*v*t)+1*1e-3*sin(2*pi/10*v.*t)+0.7*1e-3*sin(2*pi/12*v.*t) + 1e-3*sin(2*pi/30*v*t)));
fspace = @(x)(1*(2*1e-3*sin(2*pi/8 *x) + 1e-3*sin(2*pi/10*x)+0.7*1e-3*sin(2*pi/12*x) + 1e-3*sin(2*pi/30*x)));
x = 0:del_x:1600;
N = length(x);
longwave = fspace(x);
% figure1 = figure('Color',[1 1 1]);plot(x,longwave);set(gca,'Fontname','Times New Roman','fontsize',14);

%%
syms x_syms t_syms v_syms;
f_tp = diff(fspace(x_syms));
fspace_dot = @(x_syms)eval(f_tp);
wavediff = fspace_dot(x);
f_tp = diff(ftime(t_syms,v_syms) , t_syms);
f_z_v = @(t_syms,v_syms)eval(f_tp);
f_tp2 = diff(f_tp , t_syms);
f_z_acc = @(t_syms,v_syms)eval(f_tp2);

v = 100/3.6;
dt = del_x/v;                                       %%����ʱ����
t = 0:dt:(length(x)-1)*dt;
z_v = f_z_v(t,v);
% z_acc = f_z_acc(t,v) + randn(1,length(t)) * 1e-4 * 9.8;
z_acc = f_z_acc(t,v) + randn(1,length(t)) * 0 * 9.8;
% figure;plot(z_acc);hold on;plot( f_z_acc(t,v))
%%��acc��������ȫ����

%% ���ٶȼ�
z_dot =  wavediff(1)/4;%%�е�Ʈ
z_dot = 6.921463910193095e-04;
z = 0;
for i = 1:length(z_acc)
    z_dot = z_dot + z_acc(i)*dt*dt;
    z = z + z_dot;
    z_dot_save(i) = z_dot;
    z_save(i) = z;
end

% figure;plot(x,z_save*1e3);hold on;plot(x , longwave*1e3);legend 1 2
l1 = z_save-longwave;
% figure;plot(l1*1e3);

%% ���ٶȾ��������˲��������ҽ���������
%%���˳����˲���֮���ǻ���źŲ���Ӱ�죬�������Ǹ������ģ��������źŲ������䡣
z_ref_1 = shortwave_filter(z_acc*dt*dt)*5;%%30m
z_ref_2 = longwave_filter(z_acc*dt*dt,111,29,111,193)*5;%%30m
z_ref_3 = longwave_filter(z_acc*dt*dt,897,245,897,1561)*5;%%200m

figure1 = figure('Color',[1 1 1]);
% plot(z_ref_1(1:end)*5);hold on;plot( longwave(20:end) );
% plot(z_ref_3);hold on;plot( longwave(30:end) );
plot(z_ref_3(1:end));hold on;plot( longwave );
set(gca,'Fontname','Times New Roman','fontsize',16);
legend �̲��˲��� �ο���ƽ˳ �����˲���

l1_1 = z_ref_1(1:end-19)' - longwave(20:end);
figure;plot(l1_1);

%% ����fdatool�˲���֮��
% b = load('filter1.mat');
% tmpN = (length(b.Num)-1)/2;
% pmcol_70m_tmp = conv(b.Num,z_acc*dt*dt);
% pmcol_70m_tmp(1:tmpN) = [];pmcol_70m_tmp(end-tmpN+1:end) = [];
% z_dot = 3.394014695925989e-04;    %%ע�����z_dot�����������Ҫ��֤һ�ײ�ֵĳ�ֵ��û����ģ��ſ��Ա�֤���յĽ������ô֮ǰһֱ��ɢ�������ֵ�кܴ��ϵ
% zL_70m_fdatool = 0;
% for i = 1:length(pmcol_70m_tmp)
%     z_dot = z_dot + pmcol_70m_tmp(i);
%     zL_70m_fdatool = zL_70m_fdatool + z_dot;
%     z_dot_save(i) = z_dot;
%     zL_70m_fdatool_save(i) = zL_70m_fdatool;
% end
% figure1 = figure('Color',[1 1 1]);plot(zL_70m_fdatool_save+0.01);hold on;plot(longwave);
% legend fdatool��Ƶ��˲��� �������˲���;set(gca,'Fontname','Times New Roman','fontsize',16);xlabel('��� /0.25m');ylabel('�ߵ� / (32768/10 inch)');title('70m')

%% ������
% wavediff_k_1 %%w*del_t
pitch = atan(wavediff);
wy = ( pitch(2:end) - pitch(1:end-1) )/dt;  %%�����ǲ����������
wy = [0,wy ];                               %%���Ӱ��ܶ�

wy = wy + randn(1,length(wy)) * 2.424068405547680e-06;
% wy = wy + randn(1,length(wy)) * 2.424068405547680e-05;
z_dot = wavediff(1)/4;
z = 0;
for i = 1:length(wy)
    z_dot = z_dot + wy(i) * dt * del_x;
    z = z + z_dot;
    z_dot_save(i) = z_dot;
    z_save(i) = z;
end
l2 = z_save - longwave;
figure1 = figure('Color',[1 1 1]);plot(l2*1e3-0.3,'r');hold on;plot(l1*1e3,'k');legend ���ݼ�����ײ�� ���ٶȼƼ�����ײ��;
xlabel('������ /0.25m');ylabel('��� /mm')
set(gca,'Fontname','Times New Roman','fontsize',16);
%%����������˵������Ӧ������ȵģ�û��̫��Ĳ���
%%������������Ĳ����ǻ��ִ�������
%%���������Ϊ�����궨�Ĳο�����


%%
figure1 = figure('Color',[1 1 1]);
plot(([z_save]- 0.32*1e-3)*1e3);hold on;plot( longwave*1e3);
legend matlab true
%%
% figure;plot(l2*1e3);
plot_mag(longwave,'��ƽ˳Ƶ��');
% plot_mag

%% function
function plot_mag(signal_data , tit , varargin)

if (nargin == 3)
    mode = varargin{1};
    if mode == 'hold'
        hold on;
    end
else
    figure1 = figure('Color',[1 1 1]);
end

fs = 4;     %% 0.25mΪһ���������
N = length(signal_data);
% N = 1024
x = (1:N/2+1)/N*fs;
x = 1./x;
tp = abs(fftshift(fft(signal_data,N)));
tp = 20*log10(tp);
semilogx( x , tp((length(tp)/2):end)   );
xlabel('\lambda m')
ylabel('Mag dB')
set(gca,'Fontname','Times New Roman','fontsize',16);
title(tit);
grid on;
end

%% �̲��˲���
function out = shortwave_filter(in)
amcol = in;
%% 25m�����˲��ӻ���
alu = 0;elupp = 0;elup = 0;elu = 0;als = 0;alss = 0;alsss = 0;
sscal = 0.000825;
sbsci = 0.019802;
fscal = 0.1;
sbsc = 101.000;
% �����趨
Num = 768;
amcol_array = zeros(Num,1);
amcol_arraytmp = zeros(Num,1);
in1 = 533; in2 = 432;in4=382;in6=331;in7=230;
in = 539;

for i = 1:length(amcol) 
    amcol_array(in) = amcol(i);
    alu = alu + amcol_array(in1) - 3*amcol_array(in2) + 3*amcol_array(in6) - amcol_array(in7);
    elupp = alu;
    elup = elup + elupp;
    elu = elu + elup;
    emco = - amcol_array(in4);
    als = als + amcol_array(in2) - amcol_array(in6);
    alss = alss + als;
    alss = alss + sbsc*emco;
    alsss = alsss + alss;
    xtemp = (alsss*sbsci - sscal*elu)*fscal;
    yL(i,1) = xtemp;
    
    %% ��������
    in1 = mod(in1,Num)+1;
    in2 = mod(in2,Num)+1;
    in4 = mod(in4,Num)+1;
    in6 = mod(in6,Num)+1;
    in7 = mod(in7,Num)+1;
    in = mod(in,Num)+1;
    %%
    save(i,1) = alu;
    save(i,2) = elu;
    save(i,3) = als;
    save(i,4) = alsss;
end
%%
out = yL;
end
