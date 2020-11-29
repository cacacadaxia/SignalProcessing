clear all
close all;
det = 0.01;
f = 0.01:det:1e3;
lamda = 1./f;
Omega = 2*pi*f;
Cs = (1 )./( 1j*Omega);

figure1 = figure('Color',[1 1 1]);
semilogx( f , 20*log10(Cs),'LineWidth',1 );
xlabel('模拟频率 (Hz)');
ylabel('dB');
psi = 0.01:det:1e3;
fs = 100;
z_1 = exp(-j*2*pi*psi/fs);
Rz = 1./fs./( 1 - z_1 );
hold on;semilogx(psi,20*log10(Rz),'r','LineWidth',2);
set(gca,'Fontname','Times New Roman','fontsize',14);
% Rz2 = 1/fs*(1+z_1) ./ 2 ./ ( 1 - z_1 );
% semilogx(psi,20*log10(Rz2),'k','LineWidth',1);
% legend 1/s 反向差分 双线性变换

%% 

% figure;
Hz_chafen2 = 1 - 2 .* z_1 + z_1.^2;
semilogx(psi,20*log10(Hz_chafen2));


