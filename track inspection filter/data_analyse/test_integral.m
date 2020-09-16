

d = textread('gplpe.txt');
for i = 1:length(d)
    if i == 1
       x_dot = 0;
       x_dot_2 = 0
       p_x_dot = 0;
    else
        x_dot = d(i) - d(i-1);
        x_dot_2 = x_dot - p_x_dot;
        p_x_dot = x_dot;
        
        x_dot_2_s(i) = x_dot_2;
    end
end

x_dot = 0;
p_x_dot = 0;
x = 0;
p_x = 0;
for i = 4:length(d)
    
    x_dot_2 = x_dot_2_s(i);
    x_dot = x_dot + x_dot_2;
    x = x + x_dot;

    
    x_s(i) = x;
    
end

figure;
plot([zeros(1,0), x_s]);
hold on;
plot(d)

