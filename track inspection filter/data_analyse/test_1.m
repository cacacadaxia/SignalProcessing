clear;
%初始参数
L = 25;
E = 2.89e9;
sigma = 0.2;
I = 2.9;
m = 2303;
k_v = 1595;
M_v = 5750;
v = 100/3.6;

%生成参数
rho_k = 2 * k_v / (m*L);
rho_F = 2 * M_v * 9.81 / (m*L);

M = diag([1,1,1,1,1,M_v]);%质量矩阵

T = L / v;%通过时长
nt = 1000;%步数
dt = T / nt;%步长

gama = 0.25;beta = 0.5;

a0 = 1 / (gama * dt^2);
a1 = beta / (gama * dt);
a2 = 1 / (gama * dt);
a3 = 1 / (2 * gama) - 1;
a4 = beta / gama -1 ;
a5 = (beta / (2 * gama) - 1) * dt;
a6 = dt * (1 - beta);
a7 = dt * beta;

d = zeros(6,nt);%初始化广义位移，6自由度，q1-q5,z_t
v = zeros(6,nt);%速度
a = zeros(6,nt);%加速度

%newmark-beta法
for i = 2:nt %第1步已赋初值，从第2步开始到1000步
    
    t = (i - 1) * dt; %响应时间
    
    K_0 = [w_(1)^2 + rho_k * phy_(1,t) * phy_(1,t),rho_k * phy_(1,t) * phy_(2,t),...
           rho_k * phy_(1,t) * phy_(3,t),rho_k * phy_(1,t) * phy_(4,t),...
           rho_k * phy_(1,t) * phy_(5,t),-rho_k * phy_(1,t);%第一行
           
           rho_k * phy_(2,t) * phy_(1,t),w_(2)^2 + rho_k * phy_(2,t) * phy_(2,t),...
           rho_k * phy_(2,t) * phy_(3,t),rho_k * phy_(2,t) * phy_(4,t),...
           rho_k * phy_(2,t) * phy_(5,t),-rho_k * phy_(2,t);%第二行
           
           rho_k * phy_(3,t) * phy_(1,t),rho_k * phy_(3,t) * phy_(2,t),...
           w_(3)^2 + rho_k * phy_(3,t) * phy_(3,t),rho_k * phy_(3,t) * phy_(4,t),...
           rho_k * phy_(3,t) * phy_(5,t),-rho_k * phy_(3,t);%第三行
           
           rho_k * phy_(4,t) * phy_(1,t),rho_k * phy_(4,t) * phy_(2,t),...
           rho_k * phy_(4,t) * phy_(3,t),w_(4)^2 + rho_k * phy_(4,t) * phy_(4,t),...
           rho_k * phy_(4,t) * phy_(5,t),-rho_k * phy_(4,t);%第四行
           
           rho_k * phy_(5,t) * phy_(1,t),rho_k * phy_(5,t) * phy_(2,t),...
           rho_k * phy_(5,t) * phy_(3,t),rho_k * phy_(5,t) * phy_(4,t),...
           w_(5)^2 + rho_k * phy_(5,t) * phy_(5,t),-rho_k * phy_(5,t);%第五行
           
           -k_v * phy_(1,t),-k_v * phy_(2,t),-k_v * phy_(3,t),-k_v * phy_(4,t),-k_v * phy_(5,t),k_v;%第六行
           ];%初始广义刚度矩阵6*6
       
    K = K_0 + a0 * M;%用于计算的刚度矩阵
    
    F_0 = [rho_F * phy_(1,t);rho_F * phy_(2,t);rho_F * phy_(3,t)
           rho_F * phy_(4,t);rho_F * phy_(5,t);0];%广义载荷矩阵
    
    F = F_0 + M * (a0 * d(:,i-1) + a2 * v(:,i-1) + a3 * a(:,i-1));%用于计算的载荷矩阵
    
    d(:,i) = K \ F;%计算位移
    a(:,i) = a0 * (d(:,i) - d(:,i-1)) - a2 * v(:,i-1) - a3 * a(:,i-1);%加速度
    v(:,i) = v(:,i-1) + a6 * a(:,i-1) + a7 * a(:,i);%速度

end

t0=linspace(0,T,nt);

figure
plot(t0,d(1,:) * sin(1 * pi / 2) + d(2,:) * sin(2 * pi / 2) + d(3,:) * sin(3 * pi / 2) + d(4,:) * sin(4 * pi / 2) + d(5,:) * sin(5 * pi / 2))
ylabel('位移/m');
grid on;
xlabel('时间/s');
title('跨中位移响应时程');

figure
plot(t0,a(1,:) * sin(1 * pi / 2) + a(2,:) * sin(2 * pi / 2) + a(3,:) * sin(3 * pi / 2) + a(4,:) * sin(4 * pi / 2) + a(5,:) * sin(5 * pi / 2));
ylabel('加速度/m/s^2');
grid on;
xlabel('时间/s');
title('跨中加速度响应时程');

function w = w_(n)
L=25;
E=2.89e9;
I=2.9;
m=2303;
w=(n^2)*(pi^2)*((E*I)/(m*L^4))^(1/2);
end
function phy = phy_(n,t)
v=100/3.6;
L=25;
phy=sin(n*pi*v*t/L);
end

