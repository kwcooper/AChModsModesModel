% theta model



% initial learning
AEC = [1 0]; %postsynaptic input
ACA3 = [1 1 0]; %presynaptic input

W = AEC' * ACA3;

AEC = [0 1];
ACA3 = [0 1 1];


% init params
%t = 1; 

dt = 0.005; % 5 ms
thF = 8; % 8 Hz
stepsPerCycle = ceil(((1/thF)/dt));


% Encoding   (pi/2):  PLTP = 1;  PEC = 1; PCA3 = 0;
% Retrieval (3*pi/2): PLTP = -1; PEC = 0; PCA3 = 1;
PLTP = 1; % oscilates between 1 @ pi/2 & -1 @ 3*pi/2
PEC  = 1; % oscilates between 1 @ pi/2 &  0 @ 3*pi/2
PCA3 = 0; % oscilates between 1 @ pi/2 &  0 @ 3*pi/2

step = .01;

t = 1;

% define theta functions
thEC(1) = .5 * sin(t + PEC) + .5;
thCA3(1) = .5 * sin(t + PCA3) + .5;
thLTP(1) = sin(t + PLTP);

ai(:,:,1) = [0, 0];
% run sim
for t = 2:stepsPerCycle
  
thEC(t,:) = .5 * sin(t + PEC) + .5;
thCA3(t,:) = .5 * sin(t + PCA3) + .5;
thLTP(t,:) = sin(t + PLTP);  

% activation rule
ai(:,:,t) = thEC(t-1) * AEC' + thCA3(t-1) * W * ACA3';

% learning
%W_n = W + ;

end

aiR = reshape(ai, size(ai,2), size(ai,3));
figure; plot(aiR'); %plot(ai(2,:));


%% test phasic 

if 0
  % Encoding (Extrinsic)
  PLTP = 1; PEC = 1; PCA3 = 0;
  for t = 2:stepsPerCycle
    thLTP(t,:) = sin(t + PLTP);
    thEC(t,:) = .5 * sin(t + PEC) + .5;
    thCA3(t,:) = .5 * sin(t + PCA3) + .5;
  end
  
  % Retrieval (Intrinsic)
  PLTP = -1; PEC = 0; PCA3 = 1;
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

