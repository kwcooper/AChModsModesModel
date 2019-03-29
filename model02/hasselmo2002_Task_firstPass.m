function [] = hasselmo2002_Task_firstPass()

% hasselmo_2002 theta model
% runtime parameters
X = 1; 
% [0 to 1] smaller means theta has less effect on functional strength
% paper is vague about whether different values are used for the different
% layers (CA3 vs EC), for now we'll assume that it is a single fixed
% parameter.

% params
nTSteps = 25;
nCA1cells = 2;
nCA3cells = 2;
nECcells = 2;

dt = 0.005; % 5 ms
thF = 8; % 8 Hz
stepsPerCycle = ceil(((1/thF)/dt));
phaseStep = (2*pi)/stepsPerCycle; % numnber of radians to increment theta phase by for each time steps
nTrls = 5;


% allocate arrays
a_CA1 = nan(nCA1cells,nTSteps);
a_CA3 = nan(nCA3cells,nTrls); 
a_EC  = nan(nECcells, nTrls);
w_CA3 = rand(nCA1cells, nCA3cells);
w_EC  = eye(nCA1cells, nECcells); % identity matrix p.801

phase_EC = nan(nTSteps,1);   theta_EC = nan(nTSteps,1);
phase_CA3 = nan(nTSteps,1);  theta_CA3 = nan(nTSteps,1);
phase_LTP = nan(nTSteps,1);  theta_LTP = nan(nTSteps,1);

stage = 1; 

for trl = 1:nTrls
  % initialize first timestep
  % NOTE: THESE NEED TO BE SET FOR EACH PHASE OF THE TASK, TRL
 
  t = 1;
  a_CA3(:,1) = [0; 1]; % rand(nCA3cells,1) > .5;
  a_EC(:,1)  = [0; 1]; % rand(nECcells, 1) > .5;
  
  phase_EC(1)  = 0;    theta_EC(1)  = (X/2) * sin(phase_EC(1))  + (1-(X/2));
  phase_CA3(1) = pi;   theta_CA3(1) = (X/2) * sin(phase_CA3(1)) + (1-(X/2));
  phase_LTP(1) = 0;    theta_LTP(1) = sin(phase_LTP(1));
  
  a_CA1(:,1) = ((theta_EC(1) .* w_EC) * a_EC(:,stage)) + ((theta_CA3(1) .* w_CA3) * a_CA3(:,1)); % eq 2.4 p.799
  
  tempXprod(:,:,1) = (theta_LTP(1) .* a_CA1(:,1)) * a_CA3(:,stage)';
  % run one theta cycle
  for t = 2:stepsPerCycle
    
    % phasic input
    phase_EC(t)  = phase_EC(t-1)  + phaseStep;  theta_EC(t)  = (X/2) * sin(phase_EC(t))  + (1-(X/2));  % eq 2.2 p.799
    phase_CA3(t) = phase_CA3(t-1) + phaseStep;  theta_CA3(t) = (X/2) * sin(phase_CA3(t)) + (1-(X/2));  % eq 2.3 p.799
    phase_LTP(t) = phase_LTP(t-1) + phaseStep;  theta_LTP(t) = sin(phase_LTP(t));                      % eq 2.5 p.799
    
    
    %a_CA1 = w_EC .* a_EC + w_CA3 .* a_CA3; % eq 2.1
    a_CA1(:,t ) = ((theta_EC(t) .* w_EC) * a_EC(:,stage)) + ((theta_CA3(t) .* w_CA3) * a_CA3(:,stage)); % eq 2.4 p.799
    
    tempXprod(:,:,t) = (theta_LTP(t) .* a_CA1(:,t)) * a_CA3(:,stage)';
  end
  
  % compute weight updates after each theta cycle
  dw_CA3(:,:) = sum(tempXprod,3);
  
  %(repmat(theta_LTP,nCA1cells,1) .* a_CA1) * repmat(a_CA3(:,trl),1,stepsPerCycle)'); % eq. 2.6 p.801
  
  fprintf('Trial %i\n',stage);
  if any(isnan(dw_CA3)), keyboard, end
  w_CA3, dw_CA3
  w_CA3 = w_CA3 + dw_CA3;
  w_CA3
  
  if 1
    figure;
    subplot(3,1,1);
    plot(theta_EC)
    hold on; plot(theta_CA3)
    ylabel('theta_CA3');
    legend('EC','CA3'); ylim([0 1]);
    title(['trl = ', num2str(trl)]);
    
    subplot(3,1,2);
    plot(a_CA1');
    ylabel('act CA1');
    
    subplot(3,1,3);
    plot(theta_LTP)
    ylabel('theta LTP');
    pause
  end
end

% 
% 
% %% Task
% 
% %  n trials
% w_CA3() = K * a_EC() * a_CA3(); %NEEDS ATTN< ALSO DEFINE K eq 2.8 p.803
% 
% % e trials
% w_CA3() = w_CA3() + sum(); % eq. 2.9 p.803
% 
% % c trials
%  % eq. 2.10 p. 803-4
% 
% % r choice
% a_CA1() = ; % eq. 2.11 p.804
% 

%% plotting
%keyboard



figure; 
subplot(3,1,1);
plot(theta_EC)
hold on; plot(theta_CA3)
legend('EC','CA3'); ylim([0 1]);

subplot(3,1,2);
plot(a_CA1');

subplot(3,1,3);
plot(theta_LTP)
