% =========================================================================
%
%                  ����
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 11��3��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.
%        2.
%       3. 
%--------------------------------------------------------------------------


clear all;
close all;
%%
fd = fopen('accy_accx.txt','r');
A = textscan(fd,'%s%f%s%f','Delimiter',' ','CollectOutput',1);
accy = A{2};
accx = A{4};
fd = fopen('pitch_roll.txt','r');
A = textscan(fd,'%s%f%s%f','Delimiter',' ','CollectOutput',1);
wy = A{2};
wx = A{4};
fd = fopen('yaw_accz.txt','r');
A = textscan(fd,'%s%f%s%f','Delimiter',' ','CollectOutput',1);
wz = A{2};
accz = A{4};
fd = fopen('laserCameraLeft.txt','r');
A = textscan(fd,'%s%f%s%f','Delimiter',' ','CollectOutput',1);
ldpt1 = A{2};
ldpt2 = A{4};

%%
figure;plot(accx);
hold on;
plot(accy)
plot(accz)
legend a_x a_y a_z;
still_index = 7e4:7.2e4;
%%���ٶȵ�����
norm([ mean(accx(still_index)) ; mean(accy(still_index)) ; mean(accz(still_index)) - 1 ])

figure;plot(wx);hold on;plot(wy);
plot(wz);
legend roll pitch yaw;

%%
figure;plot(wx);
figure;plot(wx/100);hold on;plot(accy);grid on;legend 1 2;
%%
roll = 0;
dt = 1/500;
for i = 1:length(wx)
    roll = roll + wx(i)*dt;
    roll_save(i) = roll;
end
figure;plot(roll_save);
