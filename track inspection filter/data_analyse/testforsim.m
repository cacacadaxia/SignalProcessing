

% =========================================================================
%
%                  ����
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��16��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.
%        2.
%       3. 
%--------------------------------------------------------------------------

%% 1 ����ߵ����ݽ��д���
%% ��ȡ����
%% ��ߵ�����
% clear all;
d = textread('gplpe.txt');
d = d/129.01;

%% �Բɼ����źŽ��д���Ȼ�����˲�
%% ���ڵ�fsʼ�ն���4Hz����Ҫע�ⳤ����ʱ��Ķ�Ӧ��ϵ
lamda = 100;        %%d
omega = 2*pi/lamda;
l = length(d);
data.time = 0:0.25:((l-1)*0.25);

for i=1:length(d)
    dtp(i) = d(i) + i/1000;
end

    %% �˲�
    % b = Num;
    % sig_out = conv(data.signals.values,b);
    % figure;
    % plot( d,'LineWidth',1 );
    % hold on;
    % plot( dtp,'LineWidth',0.5 );
    % plot(  sig_out(   ((length(b)-1)/2+1):end   ) )
    % 
    % %% �����˲���֮��
    % figure;
    % fs = 4;     %% 0.25mΪһ���������
    % N = length(sig_out);
    % x = (1:N/2+1)/N*fs;
    % x = 1./x;
    % tp = abs(fftshift(fft(sig_out)));
    % tp = 20*log10(tp);
    % semilogx( x , tp((length(tp)/2):end)   );
    % xlabel('\lambda m')
    % ylabel('Mag dB')
    % set(gca,'Fontname','Times New Roman','fontsize',16);
    % %% ��ôfs����4Hz���� lamda = 1/fs;
    % % d = d-20*[sin(omega*data.time)]';
    % % figure;
    % hold on
    % fs = 4;     %% 0.25mΪһ���������
    % N = length(dtp);
    % x = (1:N/2+1)/N*fs;
    % x = 1./x;
    % tp = abs(fftshift(fft(dtp)));
    % tp = 20*log10(tp);
    % semilogx( x , tp((length(tp)/2):end)   );
    % xlabel('\lambda m')
    % ylabel('Mag dB')
    % set(gca,'Fontname','Times New Roman','fontsize',16);


%% 2 ���������������
%% ��ȡ����
%% �����������

    % d = textread('jiguang.txt');
    % d = d/129.01;
    % l = length(d);
    % data.time = 0:0.25:((l-1)*0.25);
    % data.time = 1:l;
    % data.signals.values = d;
    % data.signals.dimensions = 2;
    % 
    % sita_bt = (d(:,2)-d(:,1))/1435;
    % figure;plot(sita_bt/pi*180);


%% 3.��ȡ���е�����
%% �������е����ݽ��з���������ǿ�����matlab�д������ģ��
load_txt
sensor_data = fmctrl_data;
figure;plot(sensor_data);


%% (1)�۲����е�Ƶ��
% for i = 1:size(sensor_data,2)
%     
%     signal_data = sensor_data(:,i);
%     figure;
%     fs = 4;     %% 0.25mΪһ���������
%     N = length(signal_data);
%     x = (1:N/2+1)/N*fs;
%     x = 1./x;
%     tp = abs(fftshift(fft(signal_data)));
%     tp = 20*log10(tp);
%     semilogx( x , tp((length(tp)/2):end)   );
%     xlabel('\lambda m')
%     ylabel('Mag dB')
%     set(gca,'Fontname','Times New Roman','fontsize',16);
% end

%% �����ź�

% 1���ٶ�
tbs = sensor_data(:,end);
tbs_s = tbs/1e5;
v = 0.9./tbs_s;
distance = 0:0.25:0.25*(length(tbs)-1);
figure;plot(distance , v);
title('�ٶ� km/h')

% 2�����ٶ�
acc_z = sensor_data(:,7);

t = 0;
for i=2:length(tbs_s)
    t(i) = t(i-1) + tbs_s(i);
end
figure;plot(t./3600);
% ���0.6Сʱ


