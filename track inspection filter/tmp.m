clear
close all
global Fs
Fs = 360;
load '118m.mat'%mit���ݿ��118������
signal = val(1,100000:111600)/200;
%% ����FIR I�����20Hz���µĵ�ͨ�˲���
fp=14; fs=18;
detap = 0.01;detas = 0.01;
[M,beta] = selectFirFilterN(fp,fs,detap,detas);
N = M+1;
w = kaiser(N,beta);
hd = FIRItypeIdealpulse(fp,fs,N,'low');
h = hd.*w';
    % ��Ƶ��˲���
omega = linspace(0,pi,512);
mag = freqz(h,[1],omega);
figure
plot(omega/(2*pi)*Fs,20*log10(abs(mag)));
title('FIR��ͨ(14hz����)�˲���Ƶ����Ӧ');
xlabel('Ƶ��');
ylabel('����(dB)');
%% ����FIR I�����8Hz���ϵĸ�ͨ�˲���
fp2 = 8; fs2=4;
detap2 = 0.01; detas2 = 0.01;
[M2,beta2] = selectFirFilterN(fp2,fs2,detap2,detas2);
N2 = M2+1;
w2 = kaiser(N2,beta2);
hd2 = FIRItypeIdealpulse(fp2,fs2,N2,'high');
h2 = hd2.*w2';
    % ��Ƶ��˲���
omega = linspace(0,pi,512);
mag = freqz(h2,[1],omega);
figure
plot(omega/(2*pi)*Fs,20*log10(abs(mag)));
title('FIR��ͨ(8Hz����)�˲���Ƶ����Ӧ');
xlabel('Ƶ��');
ylabel('����(dB)');
%% �ź��˲�
sigFiltered = filter(h,[1],signal);
sigFiltered2 = filter(h2,[1],sigFiltered);
figure
subplot(3,1,1);
plot(signal,'r');
subplot(3,1,2);
hold on
plot(sigFiltered(M/2:end),'b');
subplot(3,1,3);
hold on
plot(sigFiltered2(M/2+M2/2:end),'g');
ylim([-2,2]);
title('�ź��˲�');
%% ����Ҷ�任�����˲����Ƶ��
data = FilteredSignal;
M = length(data);
N = M*2-1;
X = fft(data,N);
f = [0:M-1]*Fs/N;
figure
Xabs = abs(fftshift(X));
plot(f(1:end/2),Xabs(M:end-M/2));
title('�˲�����ź�Ƶ��');