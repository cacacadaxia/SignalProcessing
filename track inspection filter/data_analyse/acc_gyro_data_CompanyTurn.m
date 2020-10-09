
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
%  ���ܣ� 1.ת����Ϊ��׼��λ�������в���(Ŀǰ����ļ����������)
%        2.
%        3.
%
%--------------------------------------------------------------------------
close all
clear all
%% ��ȡ����
N = 20000;
filepath = 'data/0916_1337_x/';
tmp3 = textread([filepath,'fmctrl_data_1337.txt']);
if length(tmp3)>N
    tmp3 = tmp3(1:N,:);
end

%%
tmp2 = textread([filepath , 'tmp_zhongjian_1337.txt']);
if length(tmp2)>N
    tmp2 = tmp2(1:N,:);
end
% gpxbr,hfcra,lfcrp
%% ���鵥λ����
gyroll = tmp3(:,1)/246083.654;
gypitch = tmp3(:,2)/246083.654;
gyyaw = tmp3(:,6)/246083.654;
accx = tmp3(:,3)/19748;
accy = tmp3(:,4)/13852;
accz = tmp3(:,7)/13852;
% accx = accx/3276.8;
% accy = accy/3276.8;
% accz = accz/3276.8;
TBS = tmp3(:,16)/1e5;

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
%%
figure;plot(yaw_save/pi*180);
figure;plot(pitch_save/pi*180);     %%����Ȼ��������̫��
figure;plot(roll_save/pi*180);


%% ���ٶ�
for i = 1:length(accx)
    tbs = TBS(i);
    v(i) = 0.25/tbs;
end
vtp = 0;
for i = 2:length(accx)
    accx_ = accx(i) ;
    tbs = TBS(i);
    vtp(i) = vtp(i-1) + tbs * accx_;
end
figure;plot(vtp);hold on;plot(v);legend acc true

%%
