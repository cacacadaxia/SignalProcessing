



t = 25:120;
k = 2.7;
b = -0.9;
figure;
plot(t,k*t+b,'LineWidth',1);hold on;

l_i = [25,42,70,120]
l_o = [65,111,189,319];
for i = 1:4
    plot(l_i,l_o,'ro')
end

