
% =========================================================================
%
%                  ���������е��㷨
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

clear all;
close all;


%% ��ֹƵ��
wc = 0.76;%%��ͨ�˲�������ֹƵ��Ϊwc
b = wc;%%��ͨ
a = [1,wc];
[h_low,f] = freqs(b,a);
b = [1,0];%%��ͨ
[h_high,f] = freqs(b,a);
figure1 = figure('Color',[1 1 1]);
semilogx(f,20*log10(abs(h_low)));hold on;
semilogx(f,20*log10(abs(h_high)));
semilogx(f,20*log10(abs(h_low+h_high)));

set(gca,'Fontname','Times New Roman','fontsize',14);
grid on;
xlabel('\psi Hz');ylabel('Mag dB');

%% 1/s
[h,f] = freqs(1,[1,0]);
figure1 = figure('Color',[1 1 1]);
semilogx(f,20*log10(abs(h)));hold on;
set(gca,'Fontname','Times New Roman','fontsize',14);
grid on;
xlabel('\psi Hz');ylabel('Mag dB');

