function hasselmo2002_Task_V3_saddle_working_v4()
% hasselmo_2002 theta model
tic
N = 20; % 50 @ 48 seconds
saddleMat = nan(N,N);
X_CA3 = linspace(0,2*pi,N);
X_EC = linspace(0,2*pi,N);

%figure('position', [1 1 1043 788]);   

nTrls = [5]; %[1 2 3 4]
for nTrl_i = 1:length(nTrls)
  nTrls_P = nTrls(nTrl_i);
  runOne(X_EC(N/2), X_CA3(1), nTrls_P);
  
  fprintf('%i:',length(X_EC));
  for i = 1:length(X_EC)
    EC_TP = X_EC(i);
    fprintf(' %i',find(EC_TP==X_EC));
    for j = 1:length(X_CA3)
      CA3_TP = X_CA3(j);
      [saddleMat(j,i), Ms(:,j,i)] = runOne(EC_TP, CA3_TP, nTrls_P);
    end
  end
  fprintf('\n');


% figure;
% imagesc(saddleMat)
% xlabel('EC');
% ylabel('CA3');
% xlim([0 2*pi])
% ylim([0 2*pi])  
% 
% figure; 
[X,Y] = meshgrid(X_EC,X_CA3);
% surf(X,Y,saddleMat)
% title('Performance Measure')
% xlabel('Phase EC');
% ylabel('Phase CA3');
% zlabel('M');
% xlim([0 2*pi])
% ylim([0 2*pi])

figure%(99); 
nR = 4; nC = length(nTrls); 
subplotbyind(nR,nC,1,nTrl_i); surf(X,Y,squeeze(Ms(1,:,:))); xlabel('ECph'); ylabel('CA3ph'); view([45 60]); xlim([0 2*pi]); ylim([0 2*pi]);title(['trls = ' num2str(nTrls_P)]);
subplotbyind(nR,nC,2,nTrl_i); surf(X,Y,squeeze(Ms(2,:,:))); xlabel('ECph'); ylabel('CA3ph'); view([45 60]); xlim([0 2*pi]); ylim([0 2*pi]);
subplotbyind(nR,nC,3,nTrl_i); surf(X,Y,squeeze(Ms(3,:,:))); xlabel('ECph'); ylabel('CA3ph'); view([45 60]); xlim([0 2*pi]); ylim([0 2*pi]);
subplotbyind(nR,nC,4,nTrl_i); surf(X,Y,squeeze(Ms(4,:,:))); xlabel('ECph'); ylabel('CA3ph'); view([45 60]); xlim([0 2*pi]); ylim([0 2*pi]);

%figure;surf(X,Y,squeeze(Ms(4,:,:))); xlabel('ECph'); ylabel('CA3ph'); view([45 60]); xlim([0 2*pi]); ylim([0 2*pi]);
end

toc
keyboard; 
end


function [Mhat,Ms,dW] = runOne(EC_TP, CA3_TP, nTrls_P)
%todo: 
% build saddle plot
runMidTests = true; %DONT TURN THIS OFF!
%clear all;

%% runtime parameters
p.nTrls_trn = 1; % what does this mean?; THis grows this exponentally
p.nTrls_ext = 0; % what does this mean?; THis grows this exponentally
p.nTrls_rev = nTrls_P; % what does this mean?; THis grows this exponentally

p.nCA1cells = 2;
p.nCA3cells = 3;
p.nECcells = 2;

p.dt = 0.005; % 5 ms
p.thF = 8; % 8 Hz
p.stepsPerCycle = ceil(((1/p.thF)/p.dt));
p.nTSteps = p.stepsPerCycle;
p.phaseStep = (2*pi)/p.stepsPerCycle; % numnber of radians to increment theta phase by for each time steps
p.k = 1;
p.k2 = 1;
p.lrate = 0; %0.02; % learning rate, per Ehren's suggested experiment
p2 = p; 
p2.lrate = 8; % copy for the acc updates

p.thetaScale = 1; % X from Hasselmo et al (2002), p 799
% [0 to 1] smaller means theta has less effect on functional strength
% paper is vague about whether different values are used for the different
% layers (CA3 vs EC), for now we'll assume that it is a single fixed
% parameter.

% allocate arrays
a.CA1 = nan(p.nCA1cells,p.nTSteps);
a.CA3 = nan(p.nCA3cells); 
a.EC  = nan(p.nECcells);
w.CA3 = ones(p.nCA1cells, p.nCA3cells); % .* 0 
w.EC  = eye(p.nCA1cells, p.nECcells); % identity matrix p.801
tempXprod = nan(p.nCA1cells,p.nCA3cells,p.stepsPerCycle);
dwAccE.CA3 = zeros(2,3);
dwAccR.CA3 = zeros(2,3);
dwAcc.CA3 = zeros(2,3);

verbose = 0;
plt = 0;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%              Task                %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stage = 1;
if false
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%         Initial training             %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for each trial, run the theta model 
tempXprod = nan(p.nCA1cells,p.nCA3cells,p.stepsPerCycle);
for trl = 1:p.nTrls_trn
  % initialize first timestep
  % NOTE: THESE NEED TO BE SET FOR EACH PHASE OF THE TASK, TRL
  a.CA3(:,1) = [1; 1; 0;]; % rand(p.nCA3cells,1) > .5;
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
  dw.CA3(:,:) = sum(tempXprod,3)*p.dt;
  dwAcc.CA3(:,:) = dwAcc.CA3(:,:) + sum(tempXprod,3)*p.dt;

  if any(isnan(dw.CA3)), keyboard, end
  %w.CA3, dw.CA3
  w.CA3 = update_wCA3(w.CA3,dw.CA3,p);
end

if runMidTests
a.CA3(:,1) = [1; 1; 1]; % rand(p.nCA3cells,1) > .5;
a.EC(:,1) = [0; 0];
[a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,EC_TP,CA3_TP,stage);
if plt, figure(43); plotStateVariables(theta,a,'Initial Learning'); end
end
if plt, figure(42); plotStateVariables(theta,a,'Initial Learning'); end
end
ilW = p.k2;
w.CA3 = [ilW ilW 0; 0 0 0];
w.CA3_init = w.CA3;
Ms(1) = score(a,[0 1]');

   
%% Phase Two
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%         Extinction training             %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for each trial, run the theta model 
tempXprod = nan(p.nCA1cells,p.nCA3cells,p.stepsPerCycle);
for trl = 1:p.nTrls_ext
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
  dw.CA3(:,:) = sum(tempXprod,3)*p.dt;
  dwAccE.CA3(:,:) = dwAccE.CA3(:,:) + sum(tempXprod,3)*p.dt;

  
  if any(isnan(dw.CA3)), keyboard, end
  w.CA3 = update_wCA3(w.CA3,dw.CA3,p);
end

if runMidTests
a.CA3(:,1) = [1; 1; 1]; % rand(p.nCA3cells,1) > .5;
a.EC(:,1) = [0; 0];
[a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,EC_TP,CA3_TP,stage);
end
if plt, figure(54); plotStateVariables(theta,a,'Error Reversal'); end
Ms(2) = score(a,[0 1]');

%% Phase Three
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%         Reversal training             %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%keyboard;
% for each trial, run the theta model 
tempXprod = nan(p.nCA1cells,p.nCA3cells,p.stepsPerCycle);
for trl = 1:p.nTrls_rev
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
  dw.CA3(:,:) = sum(tempXprod,3)*p.dt;
  dwAccR.CA3(:,:) = dwAccR.CA3(:,:) + sum(tempXprod,3)*p.dt;

  if any(isnan(dw.CA3)), keyboard, end
  w.CA3 = update_wCA3(w.CA3,dw.CA3,p);

end
if runMidTests
  a.CA3(:,1) = [1; 1; 1]; % rand(p.nCA3cells,1) > .5;
  a.EC(:,1) = [0; 0];
  [a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,EC_TP,CA3_TP,stage);
  if plt, plotStateVariables(theta,a,'Correct Reversal'); end
end
Ms(3) = score(a,[0 1]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Test out effects of diff. learning   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply accumulated weight changes
dwAcc.CA3 = dwAccE.CA3 + dwAccR.CA3;
w.CA3 = update_wCA3(w.CA3_init,dwAccE.CA3,p2);

if runMidTests
  a.CA3(:,1) = [1; 1; 1]; % rand(p.nCA3cells,1) > .5;
  a.EC(:,1) = [0; 0];
  [a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,EC_TP,CA3_TP,stage);
  if plt, plotStateVariables(theta,a,'Correct Reversal'); end
end
Ms(2) = score(a,[0 1]');



% apply accumulated weight changes
dwAcc.CA3 = dwAccE.CA3 + dwAccR.CA3;
w.CA3 = update_wCA3(w.CA3_init,dwAccR.CA3,p2);

if runMidTests
  a.CA3(:,1) = [1; 1; 1]; % rand(p.nCA3cells,1) > .5;
  a.EC(:,1) = [0; 0];
  [a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,EC_TP,CA3_TP,stage);
  if plt, plotStateVariables(theta,a,'Correct Reversal'); end
end
Ms(3) = score(a,[0 1]');




% apply accumulated weight changes
dwAcc.CA3 = dwAccE.CA3 + dwAccR.CA3;
w.CA3 = update_wCA3(w.CA3_init,dwAcc.CA3,p2);

if runMidTests
  a.CA3(:,1) = [1; 1; 1]; % rand(p.nCA3cells,1) > .5;
  a.EC(:,1) = [0; 0];
  [a, tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,EC_TP,CA3_TP,stage);
  if plt, plotStateVariables(theta,a,'Correct Reversal'); end
end
Ms(4) = score(a,[0 1]');

Mhat = score(a,[0 1]');  
  
  
end



function [a,tempXprod, ph, theta] = runTheta(a,tempXprod,w,p,EC_TP,CA3_TP,stage)
  
  % allocate phase and theta params
  ph.EC = nan(p.nTSteps,1);   theta.EC = nan(p.nTSteps,1);
  ph.CA3 = nan(p.nTSteps,1);  theta.CA3 = nan(p.nTSteps,1);
  ph.LTP = nan(p.nTSteps,1);  theta.LTP = nan(p.nTSteps,1);

  
  ph.EC(1)  = EC_TP;   theta.EC(1)  = (p.thetaScale/2) * sin(ph.EC(1))  + (1-(p.thetaScale/2));
  ph.CA3(1) = CA3_TP;  theta.CA3(1) = (p.thetaScale/2) * sin(ph.CA3(1)) + (1-(p.thetaScale/2));
  ph.LTP(1) = 0;       theta.LTP(1) = sin(ph.LTP(1)); 
  
  % initialize CA1 activity
  a.CA1(:,1) = ((theta.EC(1) .* w.EC) * a.EC(:,stage)) + ((theta.CA3(1) .* w.CA3) * a.CA3(:,1)); % eq 2.4 p.799
  
  %if any(a.CA1<0); keyboard; end
  tempXprod(:,:,1) = (theta.LTP(1) .* a.CA1(:,1)) * a.CA3(:,stage)';
  % run one theta cycle
  for t = 2:p.stepsPerCycle
    
    % phasic input
    ph.EC(t)  = ph.EC(t-1)  + p.phaseStep;  theta.EC(t)  = (p.thetaScale/2) * sin(ph.EC(t))  + (1-(p.thetaScale/2));  % eq 2.2 p.799
    ph.CA3(t) = ph.CA3(t-1) + p.phaseStep;  theta.CA3(t) = (p.thetaScale/2) * sin(ph.CA3(t)) + (1-(p.thetaScale/2));  % eq 2.3 p.799
    ph.LTP(t) = ph.LTP(t-1) + p.phaseStep;  theta.LTP(t) = sin(ph.LTP(t));                                        % eq 2.5 p.799
    
    
    %a.CA1 = w.EC .* a.EC + w.CA3 .* a.CA3; % eq 2.1
    a.CA1(:,t) = ((theta.EC(t) .* w.EC) * a.EC(:,stage)) + ((theta.CA3(t) .* w.CA3) * a.CA3(:,stage)); % eq 2.4 p.799
    
    %a.CA1(:,t) = min(max(a.CA1(:,t),[0;0]),[1;1]); % set to cap the ca1 activity (doesn't work)
    
    
    tempXprod(:,:,t) = (theta.LTP(t) .* a.CA1(:,t)) * a.CA3(:,stage)'; %2.9
    
  end
  
end

function wCA3 = update_wCA3(wCA3,dwCA3,p)
  wCA3 = wCA3 + p.lrate.*dwCA3;
  %wCA3 = min(wCA3,p.k);
  %wCA3 = max(wCA3,p.k2);
end

function [M] = score(a,TARG)
% Scores the CA1 output
vers = 3;
switch(vers)
  case 1
    CA1_A = mean(a.CA1,2);
    Ma = (CA1_A(1) - CA1_A(2));% / (CA1_A(1) + CA1_A(2));
    Mt = TARG(1) - TARG(2);
    M = Ma/Mt;
  case 2
    M = max(TARG'*a.CA1 - (1-TARG)' * a.CA1);
  case 3
    M = TARG'*mean(a.CA1,2) - (1-TARG)' * mean(a.CA1,2);
  case 4
    M = 0;
end
% M = sum(TARG - CA1_A);
end

function plotStateVariables(theta,a,t)
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
