




%%
clear all;
close all;
n = 10 ;% number of subintervals
h = 1/n; % mesh size
x = 0:h:1; % mesh
y = Foo(x);% calc f(x)
M = MassMat1D(x); % assemble mass
b = LoadVec1D(x,@Foo) % assemble load
Pf = M\b; % solve linear system
plot(x,Pf,x,y); % plot  projection


%%
function y = Foo(x)
y=x.*sin(x)
end


function M = MassMat1D(x)
n = length(x)-1; % number of subintervals
M = zeros(n+1,n+1); % allocate mass matrix
for i = 1:n % loop over subintervals
    h = x(i+1) - x(i); % interval length
    M(i,i) = M(i,i) + h/3; % add h/3 to M(i,i)
    M(i,i+1) = M(i,i+1) + h/6;
    M(i+1,i) = M(i+1,i) + h/6;
    M(i+1,i+1) = M(i+1,i+1) + h/3;
end
end
function b = LoadVec1D(x,f)
n = length(x)-1;
b = zeros(n+1,1);
for i = 1:n
    h = x(i+1) - x(i);
    b(i) = b(i) + f(x(i))*h/2;
    b(i+1) = b(i+1) + f(x(i+1))*h/2;
end
end