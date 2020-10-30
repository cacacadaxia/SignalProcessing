
% =========================================================================
%
%                  �Ա��ظ���
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��16��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.�����˲���֮���ٻ��֣���һ�����ظ���
%        2.
%        3.
% 
%--------------------------------------------------------------------------
close all;
clear all;
filepath = 'data/0916_1337_x/';
load_txt;
size(wave_out);
N = length(fmctrl_data);
x = 0:0.25:0.25*(N-1);
x = x/1000;

%%
tmp5 = textread([filepath,'tmp2.txt']);
if length(tmp5)>N
    tmp5 = tmp5(1:N,:);
end
%%
tmp2 = textread([filepath,'tmp_zhongjian_1337.txt']);
if length(tmp2)>N
    tmp2 = tmp2(1:N,:);
end
% gpxbr,hfcra,lfcrp


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
    rou_l(i) = rou_l(i)*2;
    rou_r(i) = rou_r(i)*2;
    rou_l_dot2(i,1) = rou_l(i) - 2*rou_l(i-1) + rou_l(i-2);
    rou_r_dot2(i,1) = rou_r(i) - 2*rou_r(i-1) + rou_r(i-2);
end

% ******************step2 ���ٶȼƵ��˲� **********************************
ay = fmctrl_data(:,5);

% ---------------------- �����˲��� --------------------
for i = 1:length(ay)
    ay_Gz(i,1) = G(ay(i) , tbs(i));
end

%% ����
% sita_b = sita_b/3276.8/180*pi;
sita_b = tmp2(:,1);
for i = 3:length(sita_b)
    sita_b_dot2(i,1) = sita_b(i) - 2*sita_b(i-1) + sita_b(i-2);
end

camo = ay_Gz - G_par .* sita_b .* tbs.^2 + ht * sita_b_dot2;
camo = -camo;
% camo = floor(camo);%%����˵����ȡ������������
camo = aln(:,1);%%���������ȡ����
%%����˵����camo��Ĳ�׼����Ҫ��ay�е������
amcol = camo + rou_l_dot2;
amcor = camo - rou_r_dot2;
%% ʵ��Ŀǰ���˲���
    % FSCAL = 0.1
    % -1.036, 0.036, -0.25, 0.25, FSCAL * 2.0
    % �����趨
    FSCAL = 0.1;
    longWaveFilterBufferSize = 2048;
    mainTriWinCoef = -1.036;
    auxTriWinCoef = 0.036;
    mainRecWinCoef = -0.25;
    auxRecWinCoef = 0.25;
    scaleCoef =  FSCAL * 2.0;
    %% ���Ǵ��;��δ������ķ���ֵ������
    %% �ο�Ŀǰ����е��������ο��ĵ���������·���������ƽ˳�켼���о���
    %% Ҫ���з�װ�ԵĴ���
    %% ��������Ĵ�����ܲ��Ǻ�׼ȷ���ߺ���
    mainTriM = 160;
    auxTriM = 40;
    mainRecM = 160;
    auxRecM = 280;
    mainTriFront = 0; auxTriFront = 0; mainTriBehind = 0; auxTriBehind = 0;%%���Ǵ�
    mainRecSum = 0; auxRecSum = 0; mainTriSum = 0; auxTriSum = 0;
    MianREC = 0; AuxREC = 0; MainTRI = 0; AuxTRI = 0;
    Num = 2048;     %%ѭ�����д洢����
    data = zeros(Num,1);
    center_pos = 1000;
    dot = 0; yL = zeros(2,1);
    for i = 1:length(amcol)
        %% ��������
        data_input_index = center_pos + auxRecM;
        data(find_index(data_input_index)) = amcol(i);
        %% ���δ�
        mainRecSum = mainRecSum + data(find_index(center_pos + mainRecM )) - data(find_index(center_pos - mainRecM - 1 ));
        auxRecSum = auxRecSum + data(find_index(center_pos + auxRecM)) -  data(find_index(center_pos - auxRecM - 1 ));
        %% ���Ǵ�
        mainTriFront = mainTriFront + (data(find_index(center_pos + mainTriM-1)) - data(find_index(center_pos-1))) / mainTriM;
        mainTriBehind = mainTriBehind + ( data(find_index(center_pos - 1)) - data(find_index(center_pos - mainTriM - 1 )) ) / mainTriM;
        
        auxTriFront = auxTriFront + (data(find_index(center_pos + auxTriM-1)) - data(find_index(center_pos-1))) / auxTriM;
        auxTriBehind = auxTriBehind + ( data(find_index(center_pos - 1)) - data(find_index(center_pos - auxTriM - 1 )) ) / auxTriM;
        
        mainTriSum = mainTriSum + mainTriFront - mainTriBehind;
        auxTriSum = auxTriSum + auxTriFront - auxTriBehind;
        %% �������
        yout = data(find_index(center_pos)) + mainTriWinCoef*mainTriSum/mainTriM + auxTriWinCoef*auxTriSum/auxTriM + ...
            mainRecWinCoef*mainRecSum/(2*mainRecM+1) + auxRecWinCoef*auxRecSum/(2*auxRecM+1);
%         yout = data(center_pos) - 1.036*mainTriSum/mainTriM + 0.036*auxTriSum/auxTriM;
    %     + auxRecWinCoef*auxRecSum/(2*auxRecM+1);
        yout_save(i) = yout;
        dot = dot + yout * 2*auxTriM/(auxTriM*2+1);
        yL(i+1) = yL(i) + dot;
        %% ����
        center_pos = center_pos + 1;
        center_pos = find_index(center_pos);
    end
yL = yL/5;

%%
y1 = yL;
y1_true = wave_out(:,3);

%%
filepath = 'data/0916_1304_s/';
load_txt;
size(wave_out);
N = length(fmctrl_data);
x = 0:0.25:0.25*(N-1);
x = x/1000;

%%
tmp5 = textread([filepath,'tmp2.txt']);
if length(tmp5)>N
    tmp5 = tmp5(1:N,:);
end
%%
tmp2 = textread([filepath,'tmp_zhongjian_1337.txt']);
if length(tmp2)>N
    tmp2 = tmp2(1:N,:);
end
% gpxbr,hfcra,lfcrp


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
    rou_l(i) = rou_l(i)*2;
    rou_r(i) = rou_r(i)*2;
    rou_l_dot2(i,1) = rou_l(i) - 2*rou_l(i-1) + rou_l(i-2);
    rou_r_dot2(i,1) = rou_r(i) - 2*rou_r(i-1) + rou_r(i-2);
end

% ******************step2 ���ٶȼƵ��˲� **********************************
ay = fmctrl_data(:,5);

% ---------------------- �����˲��� --------------------
for i = 1:length(ay)
    ay_Gz(i,1) = G(ay(i) , tbs(i));
end

%% ����
% sita_b = sita_b/3276.8/180*pi;
sita_b = tmp2(:,1);
for i = 3:length(sita_b)
    sita_b_dot2(i,1) = sita_b(i) - 2*sita_b(i-1) + sita_b(i-2);
end

camo = ay_Gz - G_par .* sita_b .* tbs.^2 + ht * sita_b_dot2;
camo = -camo;
% camo = floor(camo);%%����˵����ȡ������������
camo = aln(:,1);%%���������ȡ����
%%����˵����camo��Ĳ�׼����Ҫ��ay�е������
amcol = camo + rou_l_dot2;
amcor = camo - rou_r_dot2;
%% ʵ��Ŀǰ���˲���
    % FSCAL = 0.1
    % -1.036, 0.036, -0.25, 0.25, FSCAL * 2.0
    % �����趨
    FSCAL = 0.1;
    longWaveFilterBufferSize = 2048;
    mainTriWinCoef = -1.036;
    auxTriWinCoef = 0.036;
    mainRecWinCoef = -0.25;
    auxRecWinCoef = 0.25;
    scaleCoef =  FSCAL * 2.0;
    %% ���Ǵ��;��δ������ķ���ֵ������
    %% �ο�Ŀǰ����е��������ο��ĵ���������·���������ƽ˳�켼���о���
    %% Ҫ���з�װ�ԵĴ���
    %% ��������Ĵ�����ܲ��Ǻ�׼ȷ���ߺ���
    mainTriM = 160;
    auxTriM = 40;
    mainRecM = 160;
    auxRecM = 280;
    mainTriFront = 0; auxTriFront = 0; mainTriBehind = 0; auxTriBehind = 0;%%���Ǵ�
    mainRecSum = 0; auxRecSum = 0; mainTriSum = 0; auxTriSum = 0;
    MianREC = 0; AuxREC = 0; MainTRI = 0; AuxTRI = 0;
    Num = 2048;     %%ѭ�����д洢����
    data = zeros(Num,1);
    center_pos = 1000;
    dot = 0; yL = zeros(2,1);
    for i = 1:length(amcol)
        %% ��������
        data_input_index = center_pos + auxRecM;
        data(find_index(data_input_index)) = amcol(i);
        %% ���δ�
        mainRecSum = mainRecSum + data(find_index(center_pos + mainRecM )) - data(find_index(center_pos - mainRecM - 1 ));
        auxRecSum = auxRecSum + data(find_index(center_pos + auxRecM)) -  data(find_index(center_pos - auxRecM - 1 ));
        %% ���Ǵ�
        mainTriFront = mainTriFront + (data(find_index(center_pos + mainTriM-1)) - data(find_index(center_pos-1))) / mainTriM;
        mainTriBehind = mainTriBehind + ( data(find_index(center_pos - 1)) - data(find_index(center_pos - mainTriM - 1 )) ) / mainTriM;
        
        auxTriFront = auxTriFront + (data(find_index(center_pos + auxTriM-1)) - data(find_index(center_pos-1))) / auxTriM;
        auxTriBehind = auxTriBehind + ( data(find_index(center_pos - 1)) - data(find_index(center_pos - auxTriM - 1 )) ) / auxTriM;
        
        mainTriSum = mainTriSum + mainTriFront - mainTriBehind;
        auxTriSum = auxTriSum + auxTriFront - auxTriBehind;
        %% �������
        yout = data(find_index(center_pos)) + mainTriWinCoef*mainTriSum/mainTriM + auxTriWinCoef*auxTriSum/auxTriM + ...
            mainRecWinCoef*mainRecSum/(2*mainRecM+1) + auxRecWinCoef*auxRecSum/(2*auxRecM+1);
%         yout = data(center_pos) - 1.036*mainTriSum/mainTriM + 0.036*auxTriSum/auxTriM;
    %     + auxRecWinCoef*auxRecSum/(2*auxRecM+1);
        yout_save(i) = yout;
        dot = dot + yout * 2*auxTriM/(auxTriM*2+1);
        yL(i+1) = yL(i) + dot;
        %% ����
        center_pos = center_pos + 1;
        center_pos = find_index(center_pos);
    end
yL = yL/5;

%%
y2 = yL;
y2_true = wave_out(:,3);

%%
figure;plot(y1);hold on;plot(y2);title('�����˲���֮��')
figure;plot(y1_true(1:end));hold on;plot(y2_true(3:end));title('25�����µĲ������ݶԱ�');legend 1 2


%% ����
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

function out = find_index(in)
Num = 2048;
if mod(in,Num) == 0
    out = Num;
else
    out = mod(in,Num);
end

end
