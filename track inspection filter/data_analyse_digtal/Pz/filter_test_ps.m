% =========================================================================
%
%                  验证滤波器的准确性
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月18日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.Ps这里有问题，数值太小
%        2.
%        3. 
%        4.
%        5.
% 
%--------------------------------------------------------------------------



clear all;
N = 10000;
tmp = textread('Ps3_filter_wx.txt');
if length(tmp)>N
    tmp = tmp(1:N,:);
end
% wx,tbs,hfcra,t4,...


%%
for i =1:length(tmp)
    gpro = tmp(i,1);
    tbs = tmp(i,2);
    gpro_Bz(i) = B(gpro,tbs);
    gpro_Bz(i) = gpro;
    gpro_Pz = P3z(gpro,tbs);
    out(i) = gpro_Pz;
end

figure;
plot(out(1:end));

figure;
plot(tmp(:,3))

% figure;
% plot(tmp(:,1)/100);hold on;
% plot( tmp(2:end,3) - tmp(1:end-1,3) );legend 1 2;


%% 对比
% 
% t3_tmp = tmp(:,4);
% figure;plot(t3_tmp);


%% filter

function out = B(x_k,tbs)
persistent  y;
if isempty(y)
    y = zeros(2,1);
end
y(2) = ( y(1)*2^17 + tbs*x_k )/(2^17 + tbs );

%%
y(1) = y(2);
out = y(2);

end



%%
function out = P3z(x_k,tbs_k)


persistent x y tbs;

Omega1 = 10^5/2^17;
Omega2 = 10^5/2^14;
Omega3 = 10^5/2^18;

zeta = Omega1+Omega2+Omega3;
nang = Omega2^2 + Omega2*Omega3 + Omega2*Omega1 + Omega1*Omega3 +Omega3^2;
kai = (Omega1 + Omega3)*Omega2^2 + Omega2*Omega1*Omega3 + (Omega1 + Omega2)*Omega3^2;

if isempty(x)
    x = zeros(5,1);
    y = zeros(5,1);
    tbs = zeros(2,1);
end
x(5) = x_k;
tbs(2) = tbs_k;


% T = tbs(2)/1e5;         %%这里注意
T = tbs(2);
x_dot = x(5) - x(4);
x_dot2 = x(5) - 2*x(4) + x(3);
x_dot3 = x(5) - 3*x(4) + 3*x(3) - x(4);
x_dot4 = x(5) - 4*x(4) + 6*x(3) - 4*x(2) + x(1);
N1 = kai*T^3*x(5) + nang*T^2*x_dot + zeta*T*x_dot2 + x_dot3;

y_k_1_dot = y(4) - y(3);
y_k_2_dot = y(3) - y(2);
y_k_3_dot = y(2) - y(1);
y_k_1 = y(4);       %%n-1时刻的y
N2_1 = 2^19*(1/2+1/2^5)*tbs(2)^3*y_k_1;
N2_2 = 2^37*(1/2+1/2^5+1/2^9) * tbs(2)*( tbs(1) * y_k_1_dot + tbs(2)*y_k_1 );
N2_3 = 2^51*(1/2+1/2^5) * ( tbs(1)* ( y_k_1_dot + y_k_1_dot - y_k_2_dot ) + tbs(2)*y_k_1 );
N2_4 = 2^64*( 3*( y_k_1_dot - y_k_2_dot ) + y_k_3_dot + y_k_1 );
N2 = 2^(-64) * ( N2_1 + N2_2 + N2_3 + N2_4 );

% up = 2^64*( N1/Omega1 + N2 );
up = 2^64*( N1 + N2 );
down = 2^64 + tbs(2) *2^51 *(1/2+1/2^5) + tbs(2)^2*2^51*(1/2+1/2^5+1/2^9) + tbs(2)^3*2^19*(1/2+1/2^5) + tbs(2)^4;

y(5) = up/down;
%% 更新
x(1) = x(2);
x(2) = x(3);
x(3) = x(4);
x(4) = x(5);

y(1) = y(2);
y(2) = y(3);
y(3) = y(4);
y(4) = y(5);

tbs(1) = tbs(2);
out = y(5);
end







