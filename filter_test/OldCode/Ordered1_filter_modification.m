
% =========================================================================
%
%                  һ�׼��ٶȵ�ͨ�˲���(������˲���)�뻥���˲���������̽��
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��29��
%   ���ߣ�
%--------------------------------------------------------------------------
%  ���ܣ� 1.���ȶ��ڴ��������ݵĴ������ڿ�����˲����еĴ�����������500
%        2.
%       3. 
% 10.16
%       1. �źŵĲ���Ƶ����500Hz�����Ի��ǻ��ڲ�����
%       2. �˲������������ٶȱ仯��ʵ���Ͼ�������TBS�仯
%       3. ������һ�ֱ�ʾ
% 11.9
% 1���޸�֮���
% 
%--------------------------------------------------------------------------


%% һ�����ֿ�����˲���+�����˲���
clear all;
clc;
%----------һ�����ֿ�����˲���------------
W1 = (10^5)/(2^17);
v1= 50/3.6;v2 = 120/3.6;v3 = 350/3.6;
% T1 = 0.25/v1;T2 = 0.25/v2;T3 = 0.25/v3;
T1 = 1/500;
T2 = 1/500;
T3 = 1/500;
%semilogx(f1,20*log10(abs(h1)));xlabel('Ƶ�ʣ�Hz��');ylabel('��ֵ(dB)');
figure1 = figure('Color',[1 1 1]);suptitle('һ�����ֿ����+����');
b = [ W1*T1 0];
a = [ 1+W1*T1 (-1)];
[h1 f1] = freqz(b,a,800000,500);
semilogx(v1./f1,20*log10(abs(h1)));hold on;
b = [ W1*T2 0];
a = [ 1+W1*T2 (-1)];
[h1 f1] = freqz(b,a,800000,500);
semilogx(v2./f1,20*log10(abs(h1)),'g');hold on;
b = [ W1*T3 0];
a = [ 1+W1*T3 (-1)];
[h1 f1] = freqz(b,a,800000,500);
semilogx(v3./f1,20*log10(abs(h1)),'r');hold on;
xlabel('������m��');ylabel('��ֵ(dB)');


%%
%-----------һ�����ֲ����˲���----------------
% T1 = 0.25/v1;
bc = [ 1+(W1*T1/2),  (-1)+(W1*T1/2)];
ac = [W1*T1 0];
[h1c f1c] = freqz(bc,ac,800000,500);
semilogx(v1./f1c,20*log10(abs(h1c)));hold on;

% T2 = 0.25/v2;
bc = [ 1+(W1*T2/2),  (-1)+(W1*T2/2)];
ac = [W1*T2 0];
[h1c f1c] = freqz(bc,ac,800000,500);
semilogx(v2./f1c,20*log10(abs(h1c)),'g');hold on;

bc = [ 1+(W1*T3/2),  (-1)+(W1*T3/2)];
ac = [W1*T3 0];
[h1c f1c] = freqz(bc,ac,800000,500);
semilogx(v3./f1c,20*log10(abs(h1c)),'r');hold on;


%%
%-----------------һ�����ֿ����+�����˲���-------
figure1 = figure('Color',[1 1 1]);
B1 = [2*W1*T1+(W1^2)*T1*T1 ,  (W1^2)*T1*T1-2*W1*T1];
A1 = [2*W1*T1+2*(W1^2)*T1*T1  ,  (-2)*W1*T1];
[H1 F1] = freqz(B1,A1,800000,500);
semilogx(v1./F1,20*log10(abs(H1)),'b');hold on;

B1 = [2*W1*T2+(W1^2)*T2*T2  , (W1^2)*T2*T2-2*W1*T2];
A1 = [2*W1*T2+2*(W1^2)*T2*T2  , (-2)*W1*T2];
[H1 F1] = freqz(B1,A1,800000,500);
semilogx(v2./F1,20*log10(abs(H1)),'g');hold on;

B1 = [2*W1*T3+(W1^2)*T3*T3 , (W1^2)*T3*T3-2*W1*T3];
A1 = [2*W1*T3+2*(W1^2)*T3*T3 , (-2)*W1*T3];
[H1 F1] = freqz(B1,A1,800000,500);
semilogx(v3./F1,20*log10(abs(H1)),'r');hold on;

xlabel('������m��');ylabel('��ֵ(dB)');


