function [] = hasselmo9402_working()
% hasselmo_2002 theta model


%todo: 
% add plots for ca3 and ec activity
% add ach to the model
% make a struct to hold data? 

clear all;
dbstop if error;

%% runtime parameters
p.nTrls = 10;
p.nTSteps = 25;
p.nCA1cells = 2;
p.nCA3cells = 3;
p.nECcells = 2;

p.dt = 0.005; % Default = 5 ms
p.thF = 8; % Default = 8 Hz
p.nCycles = 6; % How many cycles would we like per trial. Default = 6? 
p.stepsPerCycle = ceil(((1/p.thF)/p.dt));
p.phaseStep = (2*pi)/p.stepsPerCycle; % number of radians to increment theta phase by for each time steps
p.k = .5; % Max synaptic growth. Default = .5
p.lrate = 0.01; % learning rate, per Ehren's experiment
p.thetaScale = 1; % X from Hasselmo et al. 2002, p 799
% [0 to 1] smaller means theta has less effect on functional strength
% paper is vague about whether different values are used for the different
% layers (CA3 vs EC), for now we'll assume that it is a single fixed
% parameter. Default = 1?

% allocate arrays
a.CA1 = nan(p.nCA1cells,p.nTSteps);
a.CA3 = nan(p.nCA3cells,p.nTrls); 
a.EC  = nan(p.nECcells, p.nTrls);
w.CA3 = eye(p.nCA1cells, p.nCA3cells);
w.EC  = eye(p.nCA1cells, p.nECcells); % identity matrix p.801
tempXprod = nan(p.nCA1cells,p.nCA3cells,p.stepsPerCycle);



%% Task 

% 


stage = 1;

% for each trial, run the theta model 
for trl = 1:p.nTrls
  % initialize first timestep
  % NOTE: THESE NEED TO BE SET FOR EACH PHASE OF THE TASK, TRL
  a.CA3(:,1) = [1; 1; 0]; % rand(p.nCA3cells,1) > .5;
  a.EC(:,1)  = [1; 0]; % rand(p.nECcells, 1) > .5;
  
  % run the model
  [a, tempXprod, pha, theta] = runTheta(a,w,tempXprod,p,stage);
  
  % compute weight updates after each theta cycle
  dw.CA3(:,:) = sum(tempXprod,3);
  
  fprintf('\nTrial %i\n',trl);
  if any(isnan(dw.CA3)), keyboard, end % Just in case something goes wrong...
  w.CA3, dw.CA3
  w.CA3 = w.CA3 + p.lrate .* dw.CA3; %   w.CA3 = w.CA3 + dw.CA3;
  w.CA3 = min(w.CA3, p.k);
  %w.CA3
  
end

if 1, plotStateVariables(theta,a,pha); end

%keyboard

fprintf('Test initial learning\n');
a.CA3(:,1) = [0; 1; 0]; % rand(p.nCA3cells,1) > .5;
[a, tempXprod, pha, theta] = runTheta(a,w,tempXprod,p,stage);
if 1, plotStateVariables(theta,a,pha); ends
s

keyboard;
end 

function [a,tempXprod, pha, theta] = runTheta(a,w,tempXprod,p,stage)
  
  % allocate phase and theta params
  pha.EC = nan(p.nTSteps,1);   theta.EC = nan(p.nTSteps,1);
  pha.CA3 = nan(p.nTSteps,1);  theta.CA3 = nan(p.nTSteps,1);
  pha.LTP = nan(p.nTSteps,1);  theta.LTP = nan(p.nTSteps,1);

  
  pha.EC(1)  = 0;    theta.EC(1)  = (p.thetaScale/2) * sin(pha.EC(1))  + (1-(p.thetaScale/2));
  pha.CA3(1) = pi;   theta.CA3(1) = (p.thetaScale/2) * sin(pha.CA3(1)) + (1-(p.thetaScale/2));
  pha.LTP(1) = 0;    theta.LTP(1) = sin(pha.LTP(1));
  
  % initialize CA1 activity
  %a.CA1(:,1) = ((theta.EC(1) .* w.EC) * a.EC(:,stage)) + ((theta.CA3(1) .* w.CA3) * a.CA3(:,1)); % eq 2.4 p.799
  %p.nCycles = 10;
  %pha.AChLvls = linspace(1,0,p.nCycles*p.stepsPerCycle); % add acetylcholine (linear)
  AchPut = linspace(-10,10,p.nCycles*p.stepsPerCycle);
  pha.AChLvls = sigmoid(AchPut);
  Cr = 1; Cl = 1; % These need to be updated to control for individual dynamics in CA3 vs EC 
  updateEC3 = @(ACh,Cl,theta,w,a) (1 - ACh*Cl) * ((theta.*w) * a);
  updateCA3 = @(ACh,Cr,theta,w,a) (1 - ACh*Cr) * ((theta.*w) * a);
  
  syn.EC3(:,1) = updateEC3(pha.AChLvls(1),Cl,theta.EC(1),w.EC,a.EC(:,stage));
  syn.CA3(:,1) = updateEC3(pha.AChLvls(1),Cr,theta.CA3(1),w.CA3,a.CA3(:,stage));
%  syn.EC3(:,1) = (1 - pha.AChLvls(1)*Cl) * ((theta.EC(1) .* w.EC) * a.EC(:,stage)); % fix one minus
%  syn.CA3(:,1) = (1 - pha.AChLvls(1)*Cr) * ((theta.CA3(1) .* w.CA3) * a.CA3(:,stage));
  
  a.CA1(:,1) = syn.EC3(:,1) + syn.CA3(:,1); % eq 2.4 p.799

  tempXprod(:,:,1) = (theta.LTP(1) .* a.CA1(:,1)) * a.CA3(:,stage)';
  % run one theta cycle
  for t = 2:p.nCycles*p.stepsPerCycle
    
    % phasic input
    pha.EC(t)  = pha.EC(t-1)  + p.phaseStep;  theta.EC(t)  = (p.thetaScale/2) * sin(pha.EC(t))  + (1-(p.thetaScale/2));  % eq 2.2 p.799
    pha.CA3(t) = pha.CA3(t-1) + p.phaseStep;  theta.CA3(t) = (p.thetaScale/2) * sin(pha.CA3(t)) + (1-(p.thetaScale/2));  % eq 2.3 p.799
    pha.LTP(t) = pha.LTP(t-1) + p.phaseStep;  theta.LTP(t) = sin(pha.LTP(t));                                            % eq 2.5 p.799
    
    
    
    %a.CA1 = w.EC .* a.EC + w.CA3 .* a.CA3; % eq 2.1
    %a.CA1(:,t) = ((theta.EC(t) .* w.EC) * a.EC(:,stage)) + ((theta.CA3(t) .* w.CA3) * a.CA3(:,stage)); % eq 2.4 p.799

    syn.EC3(:,t) = (pha.AChLvls(t)*Cl) * ((theta.EC(t) .* w.EC) * a.EC(:,stage));
    syn.CA3(:,t) = (1 - pha.AChLvls(t)*Cr) * ((theta.CA3(t) .* w.CA3) * a.CA3(:,stage));
    
    a.CA1(:,t) = syn.EC3(:,t) + syn.CA3(:,t); % eq 2.4 p.799

    tempXprod(:,:,t) = (theta.LTP(t) .* a.CA1(:,t)) * a.CA3(:,stage)';
    
  end
  
  modEC = reshape(syn.EC3,2,[])';
  modCA3 = reshape(syn.CA3,2,[])';
  
  if 0
    keyboard
    
    figure; hold on;
    time = linspace(0, 1000, size(reshape(syn.EC3,2,[])',1));
    plot(time,modEC);
    plot(time,modCA3);
    xlim([0,max(time)]);
    
    figure;
    plot(time, modEC+modCA3)
    
  end
  
end

function plotStateVariables(theta,a,pha)
    figure;
    subplot(4,1,1);
    hold off; plot(theta.EC)
    hold on; plot(theta.CA3)
    ylabel('theta.CA3');
    legend('EC','CA3'); ylim([0 1]);
    title(['End of training']);
    
    subplot(4,1,2);
    plot(a.CA1');
    ylabel('act CA1');
    
    subplot(4,1,3);
    plot(theta.LTP)
    ylabel('theta LTP');
    
    subplot(4,1,4);
    plot(pha.AChLvls)
    ylabel('Ach');
end

function y = sigmoid(x)
y = 1.0 ./ (1 + exp(-x));
end
%%
