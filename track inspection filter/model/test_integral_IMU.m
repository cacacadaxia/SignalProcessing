
% =========================================================================
%
%                  轨检中的IMU数据处理的二阶差分问题讨论
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 10月15日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.生成一组数据 [m/s^2] 与 [rad/s]，IMU数据与真值
%        2.用来验证轨检的理论仿真部分
%        3.
%        4.
% 
%--------------------------------------------------------------------------

clear all;
close all;
end_t = 80;
[imu_data ,  pos , vel , time , dt]= data_gener(end_t);
PLOT = 0;
if PLOT
    figure;plot(pos(:,1:3));title('位置')
    figure;plot(pos(:,4:6));title('角度')
end

%% 加噪声
sigma_acc  = 0.00001;           % accelerometer noise (m/s^2)
sigma_gyro = 0.0*pi/180;    % gyroscope noise (rad/s)
sigma_vel = 0.0001;
%%complex add noise :randn(2,4).*[100;100000]
imu_data(1:3,:) = imu_data(1:3,:) + randn(size(imu_data(1:3,:))) * sigma_acc;
imu_data(4:6,:) = imu_data(4:6,:) + randn(size(imu_data(4:6,:))) * sigma_gyro;
vel = vel + randn(size(vel)) * sigma_acc;

%% chapter 0: 二阶差分计算
acc = vel(:,4:6);
del_t = 1/256;
s_dot = [0.010878125858769,0.02665580 , 0.01250652];
s_dot = -[0,-0.0245445353086092,-0.0122717413490623];

s = [20,5,5];
for i = 1:length(acc)
    acc_ = acc(i,:);
    tmp = acc_*del_t*del_t;
    s_dot = s_dot + tmp;
    s = s + s_dot;
    s_save(i,:) = s;
    s_dot_save(i,:) = s_dot;
end
figure1 = figure('Color',[1 1 1]);
plot( s_save - pos(:,1:3) );legend x y z
set(gca,'Fontname','Times New Roman','fontsize',14);


% acc = vel(:,4:6);
% del_t = 1/256;
% s_dot = vel(1,1:3);
% s = [20,5,5];
% for i = 1:length(acc)
%     acc_ = acc(i,:);
% %     tmp = acc_*del_t*del_t;
%     s_dot = s_dot + acc_*del_t;
%     s = s + s_dot*del_t;
%     s_save(i,:) = s;
%     s_dot_save(i,:) = s_dot;
% end
% figure;plot( s_save - pos(:,1:3) );

%%这说明了二阶差分求和与直接两次求和的一致性，在从z函数上来看两者是等价的

%% chapter 1: 中值积分
gw        = [0, 0, 9.81]'*1e3;
init_ang  = pos(1,4:6);  %%初始欧拉角
init_quat = angle2quat(init_ang(3),init_ang(2),init_ang(1));%%初始四元数
quat      = init_quat';     %%初始角度赋值
N         = length(imu_data);
Vw        = vel(1,1:3)';
Pwb       = pos(1,1:3)';
for k = 2:N
    w = (imu_data(4:6, k-1)    +    imu_data(4:6, k))/2;
    a =  imu_data(1:3, k-1);
    a1 = imu_data(1:3, k);
    Ow = [0     -w(1)   -w(2)    -w(3);...
          w(1)   0       w(3)    -w(2);...
          w(2)  -w(3)    0        w(1);...
          w(3)   w(2)   -w(1)     0  ];
    Fc = 0.5*Ow;
    F = eye(4) + Fc*dt;%%或者在这里改变积分的方法(对于矩阵的积分方法)
    quat1 = F* quat;
  
    %%
    R_S_n = quat2dcm(quat');
    R_S_n_1 = quat2dcm(quat1');
    acc_w = (R_S_n' * a + gw   +   R_S_n_1' * a1 + gw)/2;

    quat = F* quat;
    Pwb = Pwb + Vw * dt + 0.5 * dt * dt * acc_w;
    Vw = Vw + acc_w * dt;
%% data save
    quatAll(k,:) = quat';
    oula_new(k,:) = qtpoula(quat)';
    PwbSav(k,:) = Pwb';
end
%% 数据对比
real_ang = rad2deg(pos(:,4:6));          %%转换成为
figure;plot(real_ang - oula_new/pi*180);title('中值积分角度误差 /度');
figure;plot(PwbSav - pos(:,1:3));title('中值积分误差 /m');
%%直接两次积分也会出现误差？
%%直接进行中值积分，陀螺仪的角度误差会累积到位置计算中，因此这个误差比较大

%% chapter 2: kalman filter
%% 观测是速度加上误差，速度没有误差的话，自然就很准
% 观测过程
vel_measure = vel(:,1:3) + sigma_vel*randn(size(vel(:,1:3)));
% length of data
N = length(imu_data);
g = 9.81;

% set the initial state vector
x = zeros(10,1);
x(1:6) = [ pos(1,1:3)' ; vel(1,1:3)' ];
x(7:10,1) = init_quat';

% set the initial covariance
P = diag([1e-10*ones(1,6), 1e-6*ones(1,4)]);
%
x_r = zeros(10,N);
x_r(:,1) = x;

% measurement matrix
%%这个观测是不是不太妥当
H = [zeros(3), eye(3), zeros(3,4)];
R =  sigma_vel^2 * eye(3);
% 
xx = zeros(6,1);
xxsav = [];
%
%=========================================================================%
%==                             Main  Loop                               =%
%=========================================================================%
Ksave = [];
for k = 2:N
    %% compute state transition matrix F and covariance Q
    w = imu_data(4:6, k-1); % - bias_w;
    quat = x(7:10,1);
    a = imu_data(1:3, k-1); % - bias_a;
    
    % continuous state transition matrix
    %%连续状态转移矩阵
    Ow = [0     -w(1)   -w(2)    -w(3);...
          w(1)   0       w(3)    -w(2);...
          w(2)  -w(3)    0        w(1);...
          w(3)   w(2)   -w(1)     0  ];   %%q的变化速度
    Vq = compVq(quat, a);%%3,4 加速度变形之后，对于位置的速度
    
    Fc = zeros(10);
    Fc(1:3, 4:6) = eye(3);
    Fc(4:10,7:10)= [Vq; 0.5*Ow];%%10,10 连续Fc
    
    % continuous process covariance
    Gq = 0.5* [-quat(2)  -quat(3)   -quat(4); ...
                quat(1)  -quat(4)    quat(3); ...
                quat(4)   quat(1)   -quat(2); ...
               -quat(3)   quat(2)    quat(1)];       %%利用状态值更新Gq
    Qc = zeros(10);                             %%不断更新Q
    Qc(4:6, 4:6)  =  sigma_acc^2*eye(3);        %%与加速度相关的值(这里始终是相对于过程噪声来讲的)
    Qc(7:10,7:10) =  sigma_gyro^2*(Gq*Gq');     %%获得Q值的方法
    %%Gq*Gq' 4，4
    % discretilization
    F = eye(10) + Fc* dt;       %%不用管？
    Q = Qc* dt;                 %%设定Q值
    
    %% state propagation
    R_S_n = quat2dcm(quat');        %%guan to ben
    acc = R_S_n' * a + gw;  %%guan
    x(1:3) = x(1:3) + x(4:6)* dt + 0.5*acc* dt^2;
    x(4:6) = x(4:6) + acc* dt;
    
    quat = (eye(4) + 0.5*Ow* dt)*quat;%%这里的处理与加速度的不一样
    quat = quat/norm(quat);
    x(7:10) = quat;
    
    %% covariance propagation
    P = F*P*F' + Q;
    
    %% zero-velocity update
    K = (P*H')/(H*P*H'+R);
    mea = vel_measure(k,:);
    y = mea'-x(4:6);        %%这里没看懂，观测速度为0

    x = x + K*y;
    x(7:10) = x(7:10)/norm(x(7:10));

    P = (eye(10)-K*H)*P;
    P = (P+P')/2;
    x_r(:,k) = x;
    x_oula(:,k) = qtpoula(quat);
end
%% 显示数据
%%position
% figure;plot3(x_r(1,:),x_r(2,:),x_r(3,:));hold on;plot3(pos(:,1),pos(:,2),pos(:,3));
%%err
figure;plot(pos(:,1:3) - x_r(1:3,:)');title('kalman filter 位置误差 /m')
figure;plot(real_ang - x_oula'/pi*180);title('kalman filter 角度误差 /度')

if PLOT
%     figure;plot(time , PwbSav - pos(:,1:3));title('中值积分的误差对比')
end


%% fuction 
function Ang3 = qtpoula(q)
%四元数转欧拉角
x = atan2(2*(q(1)*q(2)+q(3)*q(4)),1 - 2*(q(2)^2+q(3)^2));
y = asin(2*(q(1)*q(3) - q(2)*q(4)));
z = atan2(2*(q(1)*q(4)+q(2)*q(3)),1 - 2*(q(3)^2+q(4)^2));
Ang3 = [x y z]';
end

function out = euler2Rotation(in) %%Rwb
x = in(1);
y = in(2);
z = in(3);
%欧拉角转旋转矩阵
Rx = [1      0      0;
    0 cos(x) -sin(x);
    0 sin(x) cos(x)];
Ry = [cos(y)  0 sin(y);
    0       1      0;
    -sin(y) 0 cos(y)];
Rz = [cos(z) -sin(z) 0;
    sin(z) cos(z)  0;
    0      0       1];
R = Rz*Ry*Rx;
out = R;

end

function R  =  eulerRates2bodyRates( eulerAngles)

 roll = eulerAngles(1);
 pitch = eulerAngles(2);

 cr = cos(roll);
 sr = sin(roll);
 cp = cos(pitch);
 sp = sin(pitch);

R = [1, 0, -sp;
    0, cr, sr * cp;
    0, -sr, cr * cp];
end

function out = rad2deg(in)
[h,s] = size(in);
for i = 1:h
    for j = 1:s
%         in(i,j) = mod(in(i,j) + pi,2*pi) - pi;
        in(i,j) = in(i,j)/pi*180;
        out(i,j) = mod(in(i,j) + 180, 360)-180;
    end
    
end
end
%% function 生成数据
function [imu_data ,  pos , vel ,time ,dt ] = data_gener(end_t)
dt = 1/256;
M_PI = 3.1415926;
gn = [0, 0, -9.81]*1e3;
k = 0;
time = 0:dt:end_t;
for t = 0:dt:end_t
    k = k+1;
    ellipse_x = 15;
    ellipse_y = 20;
    z = 1; 
    K1 = 10; 
    K = M_PI / 10;
    position = [ellipse_x * cos(K * t) + 5, ellipse_y * sin(K * t) + 5, z * sin(K1 * K * t) + 5];
    dp = [-K * ellipse_x * sin(K * t), K * ellipse_y * cos(K * t), z * K1 * K * cos(K1 * K * t)];
    K2 = K * K;
    ddp = [-K2 * ellipse_x * cos(K * t), -K2 * ellipse_y * sin(K * t), -z * K1 * K1 * K2 * sin(K1 * K * t)];
    %%Rotation
    k_roll = 0.1;
    k_pitch = 0.2;
    eulerAngles = [k_roll * cos(t), k_pitch * sin(t), K * t];
    eulerAnglesRates = [-k_roll * sin(t), k_pitch * cos(t), K];
    Rwb = euler2Rotation(eulerAngles);    
    imu_gyro = eulerRates2bodyRates(eulerAngles) * eulerAnglesRates';
    imu_acc = Rwb'* (ddp + gn)';

    %%
    vel(k,1:6) = [dp,ddp];
    pos(k,1:6) = [position , eulerAngles];
    imu_data(1:6,k) = [imu_acc;imu_gyro];
end
end


function Vq = compVq(q, a)

q1 = q(1);
q2 = q(2);
q3 = q(3);
q4 = q(4);

ax = a(1);
ay = a(2);
az = a(3);

v_dot_q1 = [2*ax*q1 - 2*ay*q4 + 2*az*q3;...
            2*ax*q4 + 2*ay*q1 - 2*az*q2;...
            2*ay*q2 - 2*ax*q3 + 2*az*q1];

v_dot_q2 = [2*ax*q2 + 2*ay*q3 + 2*az*q4;...
            2*ax*q3 - 2*ay*q2 - 2*az*q1;...
            2*ax*q4 + 2*ay*q1 - 2*az*q2];
 
v_dot_q3 = [2*ay*q2 - 2*ax*q3 + 2*az*q1;...
            2*ax*q2 + 2*ay*q3 + 2*az*q4;...
            2*ay*q4 - 2*ax*q1 - 2*az*q3];

v_dot_q4 = [2*az*q2 - 2*ay*q1 - 2*ax*q4;...
            2*ax*q1 - 2*ay*q4 + 2*az*q3;...
            2*ax*q2 + 2*ay*q3 + 2*az*q4];
        
Vq = [v_dot_q1, v_dot_q2, v_dot_q3, v_dot_q4];

end

