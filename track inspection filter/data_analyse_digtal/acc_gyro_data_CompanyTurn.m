
% =========================================================================
%
%                  �������ݴ���
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��24��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.ת����Ϊ��׼��λ�������в���
%        2.Ϊʲôx����ļ��ٶȾ�����ȫ��ģ�������Ᵽ��
%        3.
%   10.15
%           1. �����޸�֮��Ӧ�����ⲻ��
%
%--------------------------------------------------------------------------


%%
close all
clear all
% ��ȡ����
N = 5e5;
start_pos = 1;
filepath = 'data/0916_1337_x/';
fmctrl_data = textread([filepath,'fmctrl_data_1337.txt']);
if length(fmctrl_data)>N
    fmctrl_data = fmctrl_data(start_pos:start_pos+N-1,:);
end
tmp2 = textread([filepath , 'tmp_zhongjian_1337.txt']);
if length(tmp2)>N
    tmp2 = tmp2(start_pos : start_pos+N-1,:);
end
% gpxbr,hfcra,lfcrp

%% λ�Ƽ�
gpvlo = fmctrl_data(:, 8);
gpvro = fmctrl_data(:,10);
gpvlo = gpvlo/258.016;
gpvro = gpvro/258.016;%%mm
G = 59.5 * 25.4;
mm2o = 1/G*180/pi;
sita_bt = (gpvro-gpvlo) * mm2o;
figure1 = figure('Color',[1 1 1]);plot(sita_bt);xlabel('��� /0.25m');ylabel('\theta_{bt} /��');set(gca,'Fontname','Times New Roman','fontsize',16);

%% ���鵥λ����
gyroll = fmctrl_data(:,1)/246083.654;  %%rad/s
gypitch = fmctrl_data(:,2)/246083.654; %%rad/s
gyyaw = fmctrl_data(:,6)/246083.654;   %%rad/s
accx = fmctrl_data(:,3)/19748;%%m/s^2
accy = fmctrl_data(:,4)/13852;%%m/s^2
accz = fmctrl_data(:,7)/13852;%%m/s^2
TBS = fmctrl_data(:,16)/1e5;

figure;plot(accy);hold on;plot(accx*2);legend 1 2;
%%������һ�����⣬Ϊʲô�����ax��������ȫ�Ǵ�ģ�
%%��һ����Ҫע�⡣
%%֮����Ҫ�������ߵĽ�����жԱȡ�


%% ���ٶȼƵıȽ�

accx_body = fmctrl_data(:,13)/( 32768/2/9.8 );%%m/s^2
accy_body = fmctrl_data(:,14)/( 32768/2/9.8 );%%m/s^2
accz_body = fmctrl_data(:,15)/( 32768/2/9.8 );%%m/s^2
% figure;plot(accx);hold on; plot(accx_body);legend ����� ����
% figure;plot(accy);hold on; plot(accy_body);legend ����� ����
% figure;plot(accz);hold on; plot(accz_body);legend ����� ����
%%������������Ѿ������겨���ˣ����Լ�����ĸ�Ƶ����Ҫ��һЩ

%% ��Ƕ�
yaw = 0;pitch = 0; roll = 0;
for i = 1:length(gyyaw)
   gyyaw_ = gyyaw(i);
   gyroll_ = gyroll(i);
   gypitch_ = gypitch(i);
   tbs = TBS(i);
   yaw = yaw + tbs * gyyaw_;
   pitch = pitch + tbs*gypitch_;
   roll = roll + tbs*gyroll_;
   yaw_save(i) = yaw;
   pitch_save(i) = pitch;
   roll_save(i) = roll;
end

figure;plot(yaw_save/pi*180);title('yaw')
figure;plot(pitch_save/pi*180);title('pitch')       %%����Ȼ��������̫��
figure;plot(roll_save/pi*180);title('roll');        %%roll����һ��Ư��
%%yaw�Ƕ�ûɶ�ã�

%% ���ٶ� 
%%��ʱ���ùܣ�����������������ģ�����Ϊʲôֱ�Ӷ�acc���л�����û���õ��أ�
%%��Ϊ�����˸��ֽǶȰ��ǣ�����Ҳ���ԣ��������һֱ����x����ķ������е�
for i = 1:length(accx)
    tbs = TBS(i);
    v(i) = 0.25/tbs;
end
vtp = 0;
for i = 2:length(accx)
    accx_ = accx(i);
    tbs = TBS(i);
    vtp(i) = vtp(i-1) + tbs * accx_;
end
% figure;plot(vtp);hold on;plot(v);legend acc true

%% ����v�����˲��������˲���
b = ones(1,150)/150;
v_filter = conv(b,v);
figure;plot(v);hold on; plot(v_filter);
licheng = 0:0.25:(N-1)*0.25;
licheng = licheng/1e3;
% figure;plot(licheng,v*3.6,'r','LineWidth',1);grid on;xlabel('��� km');ylabel('�ٶ� km/h');
figure1 = figure('Color',[1 1 1]);plot(v*3.6,'r','LineWidth',1);grid on;xlabel('��� /0.25m');ylabel('�ٶ� km/h');
set(gca,'Fontname','Times New Roman','fontsize',16);

%% ����ax������ĵ�λ��m/s^2

for i = 2:length(TBS)
    accx_tbs(i) = (v_filter(i)-v_filter(i-1)) / TBS(i);
end
figure;plot(accx_tbs);




