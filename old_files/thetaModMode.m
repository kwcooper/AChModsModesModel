% simple theta mod modes model
% inspired by Hasselmo 2002

clear all;
dbstop if error;


%% initial learning

cue = 1; % set first set of associative inputs
AEC(:,cue) = [1 0]; %postsynaptic input
ACA3(:,cue) = [1 1 0]; %presynaptic input

W(:,:,1) = AEC(:,cue) * ACA3(:,cue)'; % innitial cue = 1 network


%% init params
dt = 0.005; % 5 ms
thF = 8; % 8 Hz
stepsPerCycle = ceil(((1/thF)/dt));
k = 1; % maximum synapse strength

% Phasic modulation params
PLTP = pi/2; PEC = pi/2; PCA3 = 3*pi/2;   % Encoding
%PLTP = 3*pi/2; PEC = 3*pi/2; PCA3 = pi/2; % Retrieval
%% previous activation
t = 1;

% define theta modulation functions
thEC(1) = .5 * sin(t + PEC) + .5;
thCA3(1) = .5 * sin(t + PCA3) + .5;
thLTP(1) = sin(t + PLTP);

ai(:,:,1) = [0, 0]; % set to zero for dim consistancy 

cue = 2; % set which set of associative inputs
AEC(:,cue) = [0 1];
ACA3(:,cue) = [0 1 1];

%% run model for 5 ms
for t = 2:stepsPerCycle
  
thEC(t,:) = .5 * sin(t + PEC) + .5;
thCA3(t,:) = .5 * sin(t + PCA3) + .5;
thLTP(t,:) = sin(t + PLTP);  

% calculate EC & CA3 activity
ECA(t,:) = thEC(t-1) * AEC(:,cue)';
CA3A(t,:) = thCA3(t-1) * ACA3(:,cue)';

% activation function
ai(:,:,t) = thEC(t-1) * AEC(:,cue) + thCA3(t-1) * W(:,:,t-1) * ACA3(:,cue);

% learning function
W(:,:,t) = W(:,:,t-1) + (thLTP(t-1) * ai(:,:,t)' * ACA3(:,cue)');

% adjust synapse weights to maximum
%W(W >= k) = 1;
end

% test preformance
M = AEC(:,2)' * ai(:,:,stepsPerCycle)' - AEC(:,1)' * ai(:,:,stepsPerCycle)';
%% plotting (encoding)

if 0 % plot network activity
  aiR = reshape(ai, size(ai,2), size(ai,3));
  %figure; plot(aiR'); %plot(ai(2,:));

  figure; suptitle(' ');
  subplot(3,1,1);  hold on;  plot(CA3A);  title('CA3');
  subplot(3,1,2);  hold on;  plot(ECA);   title('EC');
  subplot(3,1,3);  hold on;  plot(aiR');  title('CA1');
end

if 0 % plot weight change
  WR = reshape(W, (size(W,2) + size(W,2)), size(ai,3));
  figure; plot(WR'); title('Delta W'); xlabel('t'); ylabel('W val');
end

%% test retrieval 
if 0
  cue = 3;
  AEC(:,cue) = [0 0];
  ACA3(:,cue) = [0 0 1];
  
  for t = stepsPerCycle+1:2*stepsPerCycle
    thEC(t,:) = .5 * sin(t + PEC) + .5;
    thCA3(t,:) = .5 * sin(t + PCA3) + .5;
    thLTP(t,:) = sin(t + PLTP);
    
    ai(:,:,t) = thEC(t-1) * AEC(:,cue) + thCA3(t-1) * W(:,:,stepsPerCycle) * ACA3(:,cue);
  end
end
%% test and plot phasic modulation

if 0
  % Encoding (Extrinsic)
  PLTP = pi/2; PEC = pi/2; PCA3 = 3*pi/2;
  for t = 2:stepsPerCycle
    thLTP(t,:) = sin(t + PLTP);
    thEC(t,:) = .5 * sin(t + PEC) + .5;
    thCA3(t,:) = .5 * sin(t + PCA3) + .5;
  end
  
  % Retrieval (Intrinsic)
  PLTP = 3*pi/2; PEC = 3*pi/2; PCA3 = pi/2;
  for t = 2:stepsPerCycle
    thLTPR(t,:) = sin(t + PLTP);
    thECR(t,:) = .5 * sin(t + PEC) + .5;
    thCA3R(t,:) = .5 * sin(t + PCA3) + .5;
  end
  
  figure; suptitle(' ');
  subplot(3,1,1); hold on; plot(thLTP);  plot(thLTPR); title('LTP');
  subplot(3,1,2); hold on; plot(thEC);   plot(thECR);  title('EC');
  subplot(3,1,3); hold on; plot(thCA3);  plot(thCA3R); title('CA3');
  legend;
end

