




% 
% amcol = aln(:,3);
% yL = 0;
% xtempdot = 0;
% for i = 1:length(amcol)
%     xtempdot = amcol(i) + xtempdot;
%     xtemp = xtemp + xtempdot;
%     yL(i,1) = xtemp;
% end



doublel = camo;
intl = aln(:,1);



new = doublel - intl;
figure;plot(intl - floor(doublel));
