





% ´°º¯Êý
wc = 0.2*pi;
M = 60;
hd = ideal_lp(wc,M);
w_rect = (rectwin(M))';
h = hd .* w_rect;
fvtool(h);

function hd = ideal_lp(wc,M)
alpha = (M-1)/2;
n = [0:1:(M-1)];
m = n - alpha + eps;        % add smallest number to avoid divide by zero
hd = sin(wc*m) ./ (pi*m);   
end