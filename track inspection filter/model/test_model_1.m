



% =========================================================================
%
%                  �����㷨���ۼ���
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 10��20��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.
%        2.
%        3. 
%--------------------------------------------------------------------------

%% ��������
clear all;
close all;
del_x = 0.25;
waveLen = 300;
x = 0:del_x:waveLen;
waveMag = 40;           %%mm
longwave = waveMag*1e-3*sin(2*pi/waveLen*x);
figure;plot(x,longwave);xlabel('������� /0.25m');ylabel('��ֵ /mm');

%% ������50m��ʼ
k_1_pos = 50;
k_pos = 50 + del_x;
L = 3;%%�������������/m

wavediff_k_1 = waveMag*1e-3*2*pi/waveLen*cos(2*pi/waveLen* k_1_pos );
wavediff_k   = waveMag*1e-3*2*pi/waveLen*cos(2*pi/waveLen* k_pos );

pitch_k_1_track = atan(wavediff_k_1);
pitch_k_track =  atan(wavediff_k);
pitch_k_1 = pitch_k_1_track + 1.1e-04;  %%�������
pitch_k = pitch_k_track - 1.1e-04;

k_1_del_B = 40;
k_1_del_F = sin(pitch_k_1)*L + k_1_del_B;
k_del_B = 37;
k_del_F = sin(pitch_k)*L + k_del_B;

%% ʵ�ʵĸߵ�
 y_true = waveMag*1e-3*sin(2*pi/waveLen*k_pos) - waveMag*1e-3*sin(2*pi/waveLen*k_1_pos);
%% ����
del_fai = ( (k_del_F - k_del_B) - (k_1_del_F - k_1_del_B) )/L;
pitch_dot = pitch_k - pitch_k_1;

pitch_dot*0.25;






