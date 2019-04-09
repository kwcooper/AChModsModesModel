
% p.EC =  
% p.CA3 = 
% p.LTP = 
% X = 1;          % X = K = 1 % this should change the amplitude of the preformance measure
% 
% M = X/2 * pi * cos(p.LTP - p.EC) - K - (X/2) * pi * cos(p.LTP - p.CA3); % eq. 2.14



tic
N = 20; % 50 @ 48 seconds
saddleMat = nan(N,N);
X_CA3 = linspace(0,2*pi,N);
X_EC = linspace(0,2*pi,N);

fprintf('%i:',length(X_EC));
for i = 1:length(X_EC)
  
  fprintf(' %i',find(i==X_EC));

  for j = 1:length(X_CA3)
    
    p.EC = X_EC(i);
    p.CA3 = X_CA3(j);
    p.LTP = X_EC(i);
    X = 1; % X = K = 1 % this should change the amplitude of the preformance measure
    K = 1; % see above
    
    saddleMat(j,i) = X/2 * pi * cos(X_EC(i)) - K - (X/2) * pi * cos(X_CA3(j)); % eq. 2.14 (with phase shifts)
    %saddleMat(j,i) = X/2 * pi * cos(p.LTP - p.EC) - K - (X/2) * pi * cos(p.LTP - p.CA3); % eq. 2.14
  end
end

toc

figure;
[X,Y] = meshgrid(X_EC,X_CA3);
surf(X,Y,saddleMat); 
xlabel('ECph'); ylabel('CA3ph'); zlabel('M')
title('Performance vs phase shift')
view([45 60]); xlim([0 2*pi]); ylim([0 2*pi]);














