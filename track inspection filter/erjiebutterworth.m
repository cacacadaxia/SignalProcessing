% ����butterworth�����˲���
%������
a =10^5 / 2^14;
f  = 0.001 : 0.001 :10;  % ģ��Ƶ��Hz
omiga = f.* 2*pi;        % ģ���Ƶ�� rad/sec
Fs1 = a^2 ./((1j.*omiga).^2 + 1j.*omiga*a + a^2);

a = 10^5 / 2^13;
f  = 0.001 : 0.001 :10;  % ģ��Ƶ��Hz
omiga = f.* 2*pi;        % ģ���Ƶ�� rad/sec
Fs2 = a^2 ./((1j.*omiga).^2 + 1j.*omiga*a + a^2);

figure
plot( omiga,Fs1,'r',omiga,Fs2,'g');


% figure
% plot( omiga,20*log10(Fs1),'r',omiga,20*log10(Fs2),'g');
% xlabel('�� ');
% ylabel('dB ');
% grid on


