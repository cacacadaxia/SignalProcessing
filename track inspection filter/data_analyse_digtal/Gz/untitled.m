


clear all;
close all;
N = 10000;
aln = textread('aln.txt');
if length(aln)>N
    aln = aln(1:N,:);
end

for i =1:length(aln)
    out(i,1) = G(aln(i,1),aln(i,2));
end
figure;plot(out(30:end));hold on;
plot(aln(:,3))
legend 1 2


figure;plot(out - aln(:,3));

function out = G(x_k,tbs_k)
persistent x y1 tbs;
if isempty(x)
    x = zeros(3,1);
    y1 = zeros(2,1);
    tbs = zeros(2,1);
end

%% 更新k时刻
x(3) = x_k;
tbs(2) = tbs_k;

%% 离散滤波
x_dot = x(3) - x(2);
x_dot_2 = x(3) - 2*x(2) + x(1);
y1(2) = (x(3) + x(2)) *tbs(2)/2^15 + x_dot;
y = (y1(2)+y1(1))* (tbs(2) + tbs(1))/2^16 + x_dot_2;

%% 更新temp量
y1(1) = y1(2);
x(1) = x(2);
x(2) = x(3);
tbs(1) = tbs(2);
out = y;

end