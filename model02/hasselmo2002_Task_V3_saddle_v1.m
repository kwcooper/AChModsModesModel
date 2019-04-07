function hasselmo2002_Task_V3_saddle()
% hasselmo_2002 theta model
tic
N = 50; % 50 @ 48 seconds
saddleMat = nan(N,N);
X_CA3 = linspace(0,2*pi,N);
X_EC = linspace(0,2*pi,N);

i = 1;
j = 1;
for EC_TP = X_EC
  j = 1;
  for CA3_TP = X_CA3
    saddleMat(i,j) = runOne(EC_TP, CA3_TP);
    j = j + 1;
  end
  i = i + 1;
end

figure;
imagesc(saddleMat)
xlabel('EC');
ylabel('CA3');
xlim([0 2*pi])
ylim([0 2*pi])

figure; 
[X,Y] = meshgrid(X_EC,X_CA3);
surf(X,Y,saddleMat)
title('Performance Measure')
xlabel('Phase EC');
ylabel('Phase CA3');
zlabel('M');
xlim([0 2*pi])
ylim([0 2*pi])

toc
keyboard;
end


function [Mhat] = runOne(EC_TP, CA3_TP)
%todo: 
% build saddle plot

%clear all;

%% runtime parameters
p.nTrls = 30; % what does this mean? 
p.nTSteps = 25;
p.nCA1cells = 2;
p.nCA3cells = 3;
p.nECcells = 2;

p.dt = 0.005; % 5 ms
p.thF = 8; % 8 Hz
p.stepsPerCycle = ceil(((1/p.thF)/p.dt));
p.phaseStep = (2*pi)/p.stepsPerCycle; % numnber of radians to increment theta phase by for each time steps
p.k = .5; 
p.lrate = 0.01; % learning rate, per Ehren's suggested experiment
p.thetaScale = 1; % X from Hasselmo et al (2002), p 799
% [0 to 1] smaller means theta has less effect on functional strength
% paper is vague about whether different values are used for the different
% layers (CA3 vs EC), for now we'll assume that it is a single fixed
% parameter.

% allocate arrays
a.CA1 = nan(p.nCA1cells,p.nTSteps);
a.CA3 = nan(p.nCA3cells,p.nTrls); 
a.EC  = nan(p.nECcells, p.nTrls);
w.CA3 = zeros(p.nCA1cells, p.nCA3cells);
w.EC  = eye(p.nCA1cells, p.nECcells); % identity matrix p.801
tempXprod = nan(p.nCA1cells,p.nCA3cells,p.stepsPerCycle);

verbose = 0;
plt = 0;




%% Task 

stage = 1;


% fprintf('Probe (pre- task)\n');  
% a.CA3(:,1) = [1; 1; 1]; % rand(p.nCA3cells,1) > .5;
% a.EC(:,1) = [1; 0];
% [a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,stage); 
% if plt, plotStateVariables(theta,a,'PRE PROBE'); end


% for each trial, run the theta model 
for trl = 1:p.nTrls
  % initialize first timestep
  % NOTE: THESE NEED TO BE SET FOR EACH PHASE OF THE TASK, TRL
  a.CA3(:,1) = [1; 1; 1]; % rand(p.nCA3cells,1) > .5;
  a.EC(:,1) = [1; 0];
  
  
  % runs the model 
  [a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,EC_TP,CA3_TP,stage);
  
  if verbose
    fprintf('\nTrial %i\n',trl);
    disp('CA1 Activity');
    a.CA1
    a.CA1(1:2,end)
  end
  
  % compute weight updates after each theta cycle
  dw.CA3(:,:) = sum(tempXprod,3);
  
  if any(isnan(dw.CA3)), keyboard, end
  %w.CA3, dw.CA3
  w.CA3 = w.CA3 + p.lrate.*dw.CA3; %   w.CA3 = w.CA3 + dw.CA3;
  w.CA3 = min(w.CA3, p.k);
  %w.CA3
  
  
%   disp('SCORE:')
%   score(a)


end

if plt, plotStateVariables(theta,a,'Initial Learning'); end
%keyboard

% fprintf('Test initial learning\n');  
% a.CA3(:,1) = [1; 1; 1]; % rand(p.nCA3cells,1) > .5;
% a.EC(:,1) = [0; 1];
% [a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,stage); 
% if plt, plotStateVariables(theta,a,'TEST'); end
% disp(['SCORE:', num2str(score(a))])

% fprintf('post task probe\n');  
% a.CA3(:,1) = [1; 1; 1]; % rand(p.nCA3cells,1) > .5;
% a.EC(:,1) = [0; 0];
% [a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,stage); 
% if plt, plotStateVariables(theta,a,'RET'); end
% disp(['SCORE:', num2str(score(a))])

   
%% Phase Two

%keyboard;
% for each trial, run the theta model 
for trl = 1:p.nTrls
  a.CA3(:,1) = [1; 1; 0]; % rand(p.nCA3cells,1) > .5;
  a.EC(:,1)  = [0; 0]; % rand(p.nECcells, 1) > .5;
  
  
  % runs the model 
  [a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,EC_TP,CA3_TP,stage);
  
  if verbose
    fprintf('\nTrial %i\n',trl);
    disp('CA1 Activity');
    a.CA1
    a.CA1(1:2,end)
  end
  
  % compute weight updates after each theta cycle
  dw.CA3(:,:) = sum(tempXprod,3);
  
  
  if any(isnan(dw.CA3)), keyboard, end
  %w.CA3, dw.CA3
  w.CA3 = w.CA3 + p.lrate.*dw.CA3; %   w.CA3 = w.CA3 + dw.CA3;
  w.CA3 = min(w.CA3, p.k);
  %w.CA3
  
%   disp('SCORE:')
%   score(a)

end
  if plt, plotStateVariables(theta,a,'Error Reversal'); end 

%% Phase Three

%keyboard;
% for each trial, run the theta model 
for trl = 1:p.nTrls
  a.CA3(:,1) = [0; 1; 1]; % rand(p.nCA3cells,1) > .5;
  a.EC(:,1)  = [0; 1]; % rand(p.nECcells, 1) > .5;
  
  
  % runs the model 
  [a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,EC_TP,CA3_TP,stage);
  
  if verbose
    fprintf('\nTrial %i\n',trl);
    disp('CA1 Activity');
    a.CA1
    a.CA1(1:2,end)
  end
  
  % compute weight updates after each theta cycle
  dw.CA3(:,:) = sum(tempXprod,3);
  
  
  if any(isnan(dw.CA3)), keyboard, end
  %w.CA3, dw.CA3
  w.CA3 = w.CA3 + p.lrate.*dw.CA3; %   w.CA3 = w.CA3 + dw.CA3;
  w.CA3 = min(w.CA3, p.k);
  %w.CA3
  
%   disp('SCORE:')
%   score(a)

end
  if plt, plotStateVariables(theta,a,'Correct Reversal'); end
  

%   fprintf('Test initial learning\n');  
%   a.CA3(:,1) = [1; 1; 1]; % rand(p.nCA3cells,1) > .5;
%   a.EC(:,1) = [1; 0];
%   [a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,stage); 
%   if plt, plotStateVariables(theta,a,'TEST 2'); end
%   disp(['SCORE:', num2str(score(a))])
% 
%   fprintf('post task probe\n');  
%   a.CA3(:,1) = [1; 1; 1]; % rand(p.nCA3cells,1) > .5;
%   a.EC(:,1) = [0; 0];
%   [a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,stage); 
%   if plt, plotStateVariables(theta,a,'RET 2'); end
%   disp(['SCORE:', num2str(score(a))])

  %fprintf('post task probe\n');  
  a.CA3(:,1) = [1; 1; 1]; % rand(p.nCA3cells,1) > .5;
  a.EC(:,1) = [0; 0];
  [a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,EC_TP,CA3_TP,stage); 
  if plt, plotStateVariables(theta,a,'Choice'); end
  %disp(['SCORE:', num2str(score(a))])
  
  
  Mhat = score(a);
  
  %keyboard;
  
  
  
  
end



function [a,tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,EC_TP,CA3_TP,stage)
  
  % allocate phase and theta params
  ph.EC = nan(p.nTSteps,1);   theta.EC = nan(p.nTSteps,1);
  ph.CA3 = nan(p.nTSteps,1);  theta.CA3 = nan(p.nTSteps,1);
  ph.LTP = nan(p.nTSteps,1);  theta.LTP = nan(p.nTSteps,1);

  
  ph.EC(1)  = EC_TP;   theta.EC(1)  = (p.thetaScale/2) * sin(ph.EC(1))  + (1-(p.thetaScale/2));
  ph.CA3(1) = CA3_TP;   theta.CA3(1) = (p.thetaScale/2) * sin(ph.CA3(1)) + (1-(p.thetaScale/2));
  ph.LTP(1) = 0;    theta.LTP(1) = sin(ph.LTP(1));
  
  % initialize CA1 activity
  a.CA1(:,1) = ((theta.EC(1) .* w.EC) * a.EC(:,stage)) + ((theta.CA3(1) .* w.CA3) * a.CA3(:,1)); % eq 2.4 p.799
  
  tempXprod(:,:,1) = (theta.LTP(1) .* a.CA1(:,1)) * a.CA3(:,stage)';
  % run one theta cycle
  for t = 2:p.stepsPerCycle
    
    % phasic input
    ph.EC(t)  = ph.EC(t-1)  + p.phaseStep;  theta.EC(t)  = (p.thetaScale/2) * sin(ph.EC(t))  + (1-(p.thetaScale/2));  % eq 2.2 p.799
    ph.CA3(t) = ph.CA3(t-1) + p.phaseStep;  theta.CA3(t) = (p.thetaScale/2) * sin(ph.CA3(t)) + (1-(p.thetaScale/2));  % eq 2.3 p.799
    ph.LTP(t) = ph.LTP(t-1) + p.phaseStep;  theta.LTP(t) = sin(ph.LTP(t));                                        % eq 2.5 p.799
    
    
    %a.CA1 = w.EC .* a.EC + w.CA3 .* a.CA3; % eq 2.1
    a.CA1(:,t) = ((theta.EC(t) .* w.EC) * a.EC(:,stage)) + ((theta.CA3(t) .* w.CA3) * a.CA3(:,stage)); % eq 2.4 p.799
    
    tempXprod(:,:,t) = (theta.LTP(t) .* a.CA1(:,t)) * a.CA3(:,stage)';
    
  end
  
end

function [M] = score(a)
% Scores the CA1 output
CA1_A = mean(a.CA1,2);
M = (CA1_A(1) - CA1_A(2)) / (CA1_A(1) + CA1_A(2)); 
end

function plotStateVariables(theta,a,t)
    figure;
    subplot(3,1,1);
    hold off; plot(theta.EC)
    hold on; plot(theta.CA3)
    ylabel('theta.CA3');
    legend('EC','CA3'); ylim([0 1]);
    title(['Plot: ', t]);
    
    subplot(3,1,2);
    plot(a.CA1');
    ylabel('act CA1');
    
    subplot(3,1,3);
    plot(theta.LTP)
    ylabel('theta LTP');
end
