clear all;
%% 读取文件
fd = fopen('data.txt','r');
A = textscan(fd,'%s %f  %s  %f','CollectOutput',1);
name1 = A{1};
name2 = A{3};
data1 = A{2};
data2 = A{4};

count1 = 1;
count2 = 1;
count3 = 1;
for i=1:length(name1)   
    if name1{i} == "canpitch"
       pitch(count1) = data1(i);
       roll(count1) = data2(i);
       count1 = count1 + 1;
    elseif name1{i} == "canaccy"
        accy(count2) = data1(i);
        accx(count2) = data2(i);
        count2 = count2 + 1;
    elseif name1{i} == "canaccz"
        accz(count3) = data1(i);
        yaw(count3) = data2(i);
        count3 = count3 + 1;
    end
end

%% 分析数据
list = 284314:292314;
yaw(list) = [];


%% 标定
Fs = 500; 
omega = accy';

%% 
t0 = 1/Fs;
theta = cumsum(omega, 1)*t0;
maxNumM = 100;
L = size(theta, 1);
maxM = 2.^floor(log2(L/2));
m = logspace(log10(1), log10(maxM), maxNumM).';
m = ceil(m); % m must be an integer.
m = unique(m); % Remove duplicates.
tau = m*t0;
avar = zeros(numel(m), 1);
for i = 1:numel(m)
    mi = m(i);
    avar(i,:) = sum( (theta(1+2*mi:L) - 2*theta(1+mi:L-mi) + theta(1:L-2*mi)).^2, 1);
end
avar = avar ./ (2*tau.^2 .* (L - 2*m));
adev = sqrt(avar);
% figure
% loglog(tau, adev)
% title('Allan Deviation')
% xlabel('\tau');
% ylabel('\sigma(\tau)')
% grid on
% axis equal
% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = -0.5;
logtau = log10(tau);
logadev = log10(adev);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));
% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);
% Determine the angle random walk coefficient from the line.
logN = slope*log(1) + b;
N = 10^logN;
fprintf('角速度白噪声(连续)(或者叫角度随机游走)的标准差N为:%f\n',N);

% Plot the results.
tauN = 1;
lineN = N ./ sqrt(tau);
figure
loglog(tau, adev, tau, lineN, '--', tauN, N, 'o')
title('Allan Deviation with Angle Random Walk')
xlabel('\tau')
ylabel('\sigma(\tau)')
legend('\sigma', '\sigma_N')
text(tauN, N, 'N')
grid on
axis equal


% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = 0.5;
logtau = log10(tau);
logadev = log10(adev);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));

% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);

% Determine the rate random walk coefficient from the line.
logK = slope*log10(3) + b;
K = 10^logK;
fprintf('角速度随机游走(连续)(或者叫零偏不稳定性)的标准差K为:%f\n',K);

% Plot the results.
tauK = 3;
lineK = K .* sqrt(tau/3);
figure
loglog(tau, adev, tau, lineK, '--', tauK, K, 'o')
title('Allan Deviation with Rate Random Walk')
xlabel('\tau')
ylabel('\sigma(\tau)')
legend('\sigma', '\sigma_K')
text(tauK, K, 'K')
grid on
axis equal


% Find the index where the slope of the log-scaled Allan deviation is equal
% to the slope specified.
slope = 0;
logtau = log10(tau);
logadev = log10(adev);
dlogadev = diff(logadev) ./ diff(logtau);
[~, i] = min(abs(dlogadev - slope));

% Find the y-intercept of the line.
b = logadev(i) - slope*logtau(i);

% Determine the bias instability coefficient from the line.
scfB = sqrt(2*log(2)/pi);
logB = b - log10(scfB);
B = 10^logB;
fprintf('角速度零偏(或者叫bias、零偏)的标准差B为:%f\n',B);

% Plot the results.
tauB = tau(i);
lineB = B * scfB * ones(size(tau));
figure;
loglog(tau, adev, tau, lineB, '--', tauB, scfB*B, 'o')
title('Allan Deviation with Bias Instability')
xlabel('\tau')
ylabel('\sigma(\tau)')
legend('\sigma', '\sigma_B')
text(tauB, scfB*B, '0.664B')
grid on
axis equal

