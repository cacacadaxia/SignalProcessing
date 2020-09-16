



data = textread('test2.txt');
 figure;plot(data(:,1)/129);
 hold on
 plot(data(36:end,2)/129)
  legend 1 2
