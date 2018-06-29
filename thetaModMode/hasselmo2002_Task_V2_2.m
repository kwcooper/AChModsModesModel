function [] = hasselmo2002_Task_firstPass()
% hasselmo_2002 theta model
% runtime parameters
 
% params
nTrls = 14;
nTSteps = 25;
nCA1cells = 2;
nCA3cells = 3;
nECcells = 2;

dt = 0.005; % 5 ms
thF = 8; % 8 Hz
stepsPerCycle = ceil(((1/thF)/dt));
phaseStep = (2*pi)/stepsPerCycle; % numnber of radians to increment theta phase by for each time steps
k = .5; 
lrate = 0.01; % learning rate, per Ehren's experiment
thetaScale = 1; % X from Hasselmo et al (2002), p 799
% [0 to 1] smaller means theta has less effect on functional strength
% paper is vague about whether different values are used for the different
% layers (CA3 vs EC), for now we'll assume that it is a single fixed
% parameter.

% allocate arrays
a_CA1 = nan(nCA1cells,nTSteps);
a_CA3 = nan(nCA3cells,nTrls); 
a_EC  = nan(nECcells, nTrls);
w_CA3 = zeros(nCA1cells, nCA3cells);
w_EC  = eye(nCA1cells, nECcells); % identity matrix p.801
tempXprod = nan(nCA1cells,nCA3cells,stepsPerCycle);

stage = 1;

for trl = 1:nTrls
  % initialize first timestep
  % NOTE: THESE NEED TO BE SET FOR EACH PHASE OF THE TASK, TRL
  a_CA3(:,1) = [0; 1; 1]; % rand(nCA3cells,1) > .5;
  a_EC(:,1)  = [0; 1]; % rand(nECcells, 1) > .5;
  
  [a_CA1, tempXprod, phase_EC, phase_CA3, phase_LTP, theta_EC, theta_CA3, theta_LTP] = runTheta(a_EC, a_CA3, a_CA1, tempXprod, w_EC, w_CA3, stage, nTSteps, stepsPerCycle, phaseStep, thetaScale);
  
  % compute weight updates after each theta cycle
  dw_CA3(:,:) = sum(tempXprod,3);
  
  fprintf('\nTrial %i\n',trl);
  if any(isnan(dw_CA3)), keyboard, end
  w_CA3, dw_CA3
  w_CA3 = w_CA3 + lrate.*dw_CA3; %   w_CA3 = w_CA3 + dw_CA3;
  w_CA3 = min(w_CA3, k);  
  w_CA3
  
end

if 1, plotStateVariables(theta_EC, theta_CA3,a_CA1,theta_LTP); end

  keyboard

  fprintf('Test initial learning\n');  
  a_CA3(:,1) = [0; 1; 0]; % rand(nCA3cells,1) > .5;
  [a_CA1, tempXprod, phase_EC, phase_CA3, phase_LTP, theta_EC, theta_CA3, theta_LTP] = runTheta(a_EC, a_CA3, a_CA1,tempXprod, w_EC, w_CA3, stage, nTSteps, stepsPerCycle, phaseStep, thetaScale); 
  if 1, plotStateVariables(theta_EC, theta_CA3,a_CA1,theta_LFP); end
  
end

function [a_CA1,tempXprod, phase_EC, phase_CA3, phase_LTP, theta_EC, theta_CA3, theta_LTP] = runTheta(a_EC, a_CA3, a_CA1,tempXprod, w_EC, w_CA3, stage, nTSteps, stepsPerCycle, phaseStep, thetaScale)

  phase_EC = nan(nTSteps,1);   theta_EC = nan(nTSteps,1);
  phase_CA3 = nan(nTSteps,1);  theta_CA3 = nan(nTSteps,1);
  phase_LTP = nan(nTSteps,1);  theta_LTP = nan(nTSteps,1);


  phase_EC(1)  = 0;    theta_EC(1)  = (thetaScale/2) * sin(phase_EC(1))  + (1-(thetaScale/2));
  phase_CA3(1) = pi;   theta_CA3(1) = (thetaScale/2) * sin(phase_CA3(1)) + (1-(thetaScale/2));
  phase_LTP(1) = 0;    theta_LTP(1) = sin(phase_LTP(1));
  
  a_CA1(:,1) = ((theta_EC(1) .* w_EC) * a_EC(:,stage)) + ((theta_CA3(1) .* w_CA3) * a_CA3(:,1)); % eq 2.4 p.799
  
  tempXprod(:,:,1) = (theta_LTP(1) .* a_CA1(:,1)) * a_CA3(:,stage)';
  % run one theta cycle
  for t = 2:stepsPerCycle
    
    % phasic input
    phase_EC(t)  = phase_EC(t-1)  + phaseStep;  theta_EC(t)  = (thetaScale/2) * sin(phase_EC(t))  + (1-(thetaScale/2));  % eq 2.2 p.799
    phase_CA3(t) = phase_CA3(t-1) + phaseStep;  theta_CA3(t) = (thetaScale/2) * sin(phase_CA3(t)) + (1-(thetaScale/2));  % eq 2.3 p.799
    phase_LTP(t) = phase_LTP(t-1) + phaseStep;  theta_LTP(t) = sin(phase_LTP(t));                      % eq 2.5 p.799
    
    
    %a_CA1 = w_EC .* a_EC + w_CA3 .* a_CA3; % eq 2.1
    a_CA1(:,t) = ((theta_EC(t) .* w_EC) * a_EC(:,stage)) + ((theta_CA3(t) .* w_CA3) * a_CA3(:,stage)); % eq 2.4 p.799
    
    tempXprod(:,:,t) = (theta_LTP(t) .* a_CA1(:,t)) * a_CA3(:,stage)';
  end
  
end

function plotStateVariables(theta_EC, theta_CA3,a_CA1,theta_LFP)
    figure;
    subplot(3,1,1);
    hold off; plot(theta_EC)
    hold on; plot(theta_CA3)
    ylabel('theta_CA3');
    legend('EC','CA3'); ylim([0 1]);
    title(['End of training']);
    
    subplot(3,1,2);
    plot(a_CA1');
    ylabel('act CA1');
    
    subplot(3,1,3);
    plot(theta_LTP)
    ylabel('theta LTP');
end
