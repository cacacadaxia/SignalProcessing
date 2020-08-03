function [L, M, N] = winlen(K)

k21 = 0.59;
k31 = 1;
k41 = 1.7368;

% a1 = 1.0;
% a2 = 1.0;

% M1 = 319;%188;%112;%65;
L = round(K*k21);
M = round(K*k31);
N = round(K*k41);