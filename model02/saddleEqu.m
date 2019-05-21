
% p.EC =  
% p.CA3 = 
% p.LTP = 
% X = 1;          % X = K = 1 % this should change the amplitude of the preformance measure
% 
% M = X/2 * pi * cos(p.LTP - p.EC) - K - (X/2) * pi * cos(p.LTP - p.CA3); % eq. 2.14

N = 20; % 50 @ 48 seconds
saddleMats = nan(N,N,3);
iter = [1,2,3];
for part = iter
    tic

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

        switch(part)
            case 1
                saddleMats(j,i,part) = X/2 * pi * cos(X_EC(i));
                
            case 2 
                saddleMats(j,i,part) = - K - (X/2) * pi * cos(X_CA3(j));
                %saddleMat(j,i) = - K - (X/2) * pi * cos(X_CA3(j));

            case 3
                saddleMats(j,i,part) = X/2 * pi * cos(X_EC(i)) - K - (X/2) * pi * cos(X_CA3(j)); % eq. 2.14 (with phase shifts)
                 %saddleMat(j,i) = X/2 * pi * cos(p.LTP - p.EC) - K - (X/2) * pi * cos(p.LTP - p.CA3); % eq. 2.14
                
        end 

      end
    end

    toc

end

figure('position', [1 1 875 245]);
names = ["EC Part", "CA3 part", "Whole F(x)"];
for sPi = iter
    subplot(1,length(iter),sPi);
    [X,Y] = meshgrid(rad2deg(X_EC),rad2deg(X_CA3));
    surf(X,Y,saddleMats(:,:,sPi)); 
    xlabel('\Phi EC'); ylabel('\Phi CA3'); zlabel('M')
    title(names(sPi)); 
    view([45 60]); 
    xlim([0 360]); ylim([0 360]); zlim([-5 3]);
    xticks([0 180 360]); yticks([0 180 360]); zticks([-3 0 3]);
    %xlim([0 2*pi]); ylim([0 2*pi]);
    colormap('winter')
end













