
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
%        2.�����˲����Ĳ��֣������жԱȡ���Ҫ�ǹ��򣬽���ܺ�
%        3. �޸Ļ��ֵķ���������ʵû�У����ǰѳ����е��˲�����д��
% 
% 
%--------------------------------------------------------------------------
close all;
clear all;
load_txt;
size(wave_out);
N = length(fmctrl_data);
x = 0:0.25:0.25*(N-1);
x = x/1000;
%%
tmp5 = textread('data/0916_1337/tmp2.txt');
if length(tmp5)>N
    tmp5 = tmp5(1:N,:);
end
%%
tmp2 = textread('data/0916_1337/tmp_zhongjian_1337.txt');
if length(tmp2)>N
    tmp2 = tmp2(1:N,:);
end
% gpxbr,hfcra,lfcrp

%% ���ĶԱ�

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
    rou_l(i) = rou_l(i)*2;
    rou_r(i) = rou_r(i)*2;
    rou_l_dot2(i,1) = rou_l(i) - 2*rou_l(i-1) + rou_l(i-2);
    rou_r_dot2(i,1) = rou_r(i) - 2*rou_r(i-1) + rou_r(i-2);
end

% ******************step2 ���ٶȼƵ��˲� **********************************
ay = fmctrl_data(:,5);
gpan = fmctrl_data(:,5);


gpan = gpan/(1359.5241/0.01/9.83);
figure;plot(gpan);
x = 0;x_dot = 0;
for i = 1:length(gpan)
    x_dot = x_dot + gpan(i)*tbs(i);
    x_dot_save(i) = x_dot;
    x = x + x_dot*tbs(i);
    x_save(i) = x;
end


% ---------------------- �����˲��� --------------------
for i = 1:length(ay)
    ay_Gz(i,1) = G(ay(i) , tbs(i));
    %% ��Ȼ��������һ��������ǣ�Gz˳���*T^2�Ļ
end

%% Ƶ�׹۲�
% plot_mag(ay,'�˲�ǰ')
% plot_mag(ay_Fz,'�˲���')
% 
% figure;plot(aln(:,2)-ay_Gz);


%% ����
% sita_b = sita_b/3276.8/180*pi;
sita_b = tmp2(:,1);
for i = 3:length(sita_b)
    sita_b_dot2(i,1) = sita_b(i) - 2*sita_b(i-1) + sita_b(i-2);
end

camo = ay_Gz - G_par .* sita_b .* tbs.^2 + ht * sita_b_dot2;
%% ����������޸�
gyroYaw = fmctrl_data(:,6);

gyroYaw_p = 0;
for i = 1:length(gyroYaw)
    gyroYaw_ = gyroYaw(i);
    tbs_ = tbs(i);
    omiga1 = 0.76;
    gyawTemp1 = gyroYaw_ - gyroYaw_p;
    gyawTemp2 = (gyroYaw_ + gyroYaw_p)*tbs_/2^18;
    
    gpYaw = (gyawTemp1 + gyawTemp2) / omiga1;
    %%
    gyroYaw_p = gyroYaw_;

    yaw_Rz(i,1) = gpYaw;
    
end
sampleDistance = 0.25;
yawParameter = 2.0970;      %% 135751/9.8*6.0135*6.0135/4294.97*pi/180
gpyawReviseTemp1 = yaw_Rz*yawParameter*sampleDistance;
gpyawRevise = gpyawReviseTemp1 + 0;

camo = - gpyawRevise + ht * sita_b_dot2;
camo = -camo;
camo = zeros(size(camo));%%���������ȫ�����0
% figure;plot(camo );hold on;plot(aln(:,1));legend matlab gj

%% ֻ�ǽ��������򵥵���ֲ֮�󣬷����������Ǵ�������ģ�������ʱ�ȷ�һ��
%%
% camo = floor(camo);%%����˵����ȡ������������
% camo = aln(:,1);%%���������ȡ����
%%����˵����camo��Ĳ�׼����Ҫ��ay�е������
amcol = camo + rou_l_dot2;
amcor = camo - rou_r_dot2;
%%
marm = aln(:,5);
g = aln(:,6);

figure;plot(amcol,'LineWidth',1);hold on;plot(aln(:,3));title('amcol֮��ĶԱ�');legend matlab gj
% figure;plot(amcol-aln(:,3));title('amcol֮��ĶԱ�');

 %%
alu = 0;elupp = 0;elup = 0;elu = 0;als = 0;alss = 0;alsss = 0;
sscal = 0.000825;
sbsci = 0.019802;
fscal = 0.1;
sbsc = 101.000;
%% ��Ȼ�������ķ���ȥʵ���ˣ����ǲ�û�л����ͬ�Ľ��
%% ���鲻֪����ɶ�þ�����
%% ����Ļ�����ʲô�����أ�
%% �иĽ��Ŀռ�
Num = 768;
amcol_array = zeros(Num,1);
amcol_arraytmp = zeros(Num,1);
in1 = 533; in2 = 432;in4=382;in6=331;in7=230;%%��Щ�������ǹ̶��ģ����ܱ�
in = 539;%%����ƫ��ֵ

%%
for i = 1:length(amcol)           %%�򵥻��֣��϶��ǲ��Ե�
    %%
   amcol_array(in) = amcol(i);
    %%�ɴ˿ɼ���������������⣬���Բ�������
    %%-----------------------------
    
    alu = alu + amcol_array(in1) - 3*amcol_array(in2) + 3*amcol_array(in6) - amcol_array(in7);
    elupp = alu;
    elup = elup + elupp;
    elu = elu + elup;
    emco = - amcol_array(in4);

    als = als + amcol_array(in2) - amcol_array(in6);
    
    alss = alss + als + sbsc*emco;
    alsss = alsss + alss;
    xtemp = (sbsci*alsss - sscal*elu)*fscal;
    
    alsss_save(i) = alsss;
    elu_save(i) = elu;
    alss_save(i) = alss;
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
figure;plot(yL,'LineWidth',1);hold on;plot(aln(:,4));legend matlab gj
figure;plot((yL - aln(:,4))/103);%%������ȫһ��
title('�����Ľ��(��λ��mm)');
%%
figure;plot(alsss_save);hold on;plot(elu_save)
%%  
    % guixiang_l = wave_out(:,3);
    % figure;plot(guixiang_l(268:end),'LineWidth',1);
    % hold on;plot(aln(:,4));
    % title('�����˲��Ľ��');

%% ��������Ļ��ֲ���Ҫ����ʱ������
    % for i=1:length(amcol)
    %     yL_2(i) = integrational(amcol(i));
    % end
    
%% �����˲���










%% ǰ���˲���
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



%% �ڳ����˲�֮����л��֣�
function out = integrational(auat)
persistent  euas euasp  auatp eusa;
if isempty(euasp)
    euasp =  0;
    auatp = 0;
    eusa = 0;
    euas = 0;
end
gyroabase = 57;

%% ����
euas = euas + (auat + auatp) * 0.5;
eusa = (2 * gyroabase / (2 * gyroabase + 1) * eusa) + (euas + euasp) * 0.5;

%% ����
out = eusa;
auatp = auat;
euasp = euas;

end





% %% �����������
% 
% for i = 1:length(amcol)           %%�򵥻��֣��϶��ǲ��Ե�
%     %%
%    amcol_array(in) = amcol(i);
%     %%�ɴ˿ɼ���������������⣬���Բ�������
%     %%-----------------------------
%     
%     alu = alu + amcol_array(in1) - 3*amcol_array(in2) + 3*amcol_array(in6) - amcol_array(in7);
%     elupp = alu;
%     elup = elup + elupp;
%     elu = elu + elup;
%     emco = - amcol_array(in4);
% 
%     als = als + amcol_array(in2) - amcol_array(in6);
%     
%     alss = alss + als + sbsc*emco;
%     alsss = alsss + alss;
%     xtemp = (alsss*sbsci - sscal*elu)*fscal;
%     yL(i,1) = xtemp;
%     
%     %% ��������
%     in1 = mod(in1,Num)+1;
%     in2 = mod(in2,Num)+1;
%     in4 = mod(in4,Num)+1;
%     in6 = mod(in6,Num)+1;
%     in7 = mod(in7,Num)+1;
%     in = mod(in,Num)+1;
%     %%
%     save(i,1) = alu;
%     save(i,2) = elu;
%     save(i,3) = als;
%     save(i,4) = alsss;
% end

