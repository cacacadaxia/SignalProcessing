mag = abs(Acc_Trap_42);

%% find index
index = find_707(mag);
ind = find(index==1);
horiz = mag(ind(1));
horiz_x = lamda(ind(1));
verti = 25;

%%  plot
% figure;
semilogx((lamda), ((mag)),'LineWidth',1);
hold on;
semilogx([1,1000] , [horiz,horiz] ,'--','LineWidth',1)
semilogx( [horiz_x,horiz_x],[0.3,0.9],'--','LineWidth',1)
semilogx( [verti,verti],[0.3,0.9],'--','LineWidth',1)


%% function
function index = find_707(mag)
tol = 0.00;
index = zeros(1,length(mag));
while(sum(index)==0)
    index1 = mag <0.7079+tol;
    index2 = mag >0.7079-tol;
    index = index1.*index2;
    tol = tol+0.001;
end
end



