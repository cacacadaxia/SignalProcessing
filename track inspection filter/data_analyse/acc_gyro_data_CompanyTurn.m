
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
%        2.Ϊʲôx����ļ��ٶȾ�����ȫ��ģ�
%        3.
%
%--------------------------------------------------------------------------
close all
clear all

%% ��ȡ����
N = 5e4;
start_pos = 14;
filepath = 'data/0916_1337_x/';
tmp3 = textread([filepath,'fmctrl_data_1337.txt']);
if length(tmp3)>N
    tmp3 = tmp3(start_pos:start_pos+N-1,:);
end

%%
tmp2 = textread([filepath , 'tmp_zhongjian_1337.txt']);
if length(tmp2)>N
    tmp2 = tmp2(start_pos : start_pos+N-1,:);
end
% gpxbr,hfcra,lfcrp

%% ���鵥λ����
gyroll = tmp3(:,1)/246083.654;  %%rad/s
gypitch = tmp3(:,2)/246083.654; %%rad/s
gyyaw = tmp3(:,6)/246083.654;   %%rad/s
accx = tmp3(:,3)/19748;%%m/s^2
accy = tmp3(:,4)/13852;%%m/s^2
accz = tmp3(:,7)/13852;%%m/s^2
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
figure;plot(yaw_save/pi*180);title('yaw')
figure;plot(pitch_save/pi*180);title('pitch')       %%����Ȼ��������̫��
figure;plot(roll_save/pi*180);title('roll');        %%roll����һ��Ư��


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
figure;plot(licheng,v*3.6,'r','LineWidth',1);grid on;xlabel('��� km');ylabel('�ٶ� km/h');
set(gca,'Fontname','Times New Roman','fontsize',16);



%% ����ax������ĵ�λ��m/s^2

for i = 2:length(TBS)
    accx_tbs(i) = (v_filter(i)-v_filter(i-1)) / TBS(i);
end
% figure;plot(accx_tbs);


%%

for i = 2:length(TBS)
    gyro_tmp(i) = 0.5*yaw_save(i)*(TBS(i) - TBS(i-1))/( TBS(i) + TBS(i-1) );
end
