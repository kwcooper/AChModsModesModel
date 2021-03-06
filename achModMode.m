
% Hasselmo ACh model

%% Params and allocation

numCells = 3;
tS = 4;

% post/pre syn input    % activation
Ai = nan(numCells,tS); ai = nan(numCells,tS+1);
Aj = nan(numCells,tS); 
W = nan(numCells,numCells, tS); % synapses
c = [1 1 1 1 1]; % ACh input

%% Run model

% initialize mats 
Ai(:,1) = [1 0 1]'; 
Ai(:,2) = [0 1 0]'; 
Ai(:,3) = [0 1 0]'; 
Ai(:,4) = [0 0 0]';

Aj(:,1) = [1 1 0]';
Aj(:,2) = [0 1 1]';
Aj(:,3) = [0 1 1]';
Aj(:,4) = [0 0 1]';

ai(:,1) = nan(numCells,1);
W(:,:,1) = zeros(numCells,numCells);

for t = 2:tS
  ai(:,t) = Ai(:,t-1) + c(t-1) * W(:,:,t-1) * Aj(:,t-1); % activity computation
  W(:,:,t) = W(:,:,t-1) + ai(:,t) * Aj(:,t-1)'; % Learning rules
end  

keyboard;
%%
% activity computation
%ai(2) = Ai(1) + c(1) * W * Aj(1);
% example
%ai = [1 0 1]' + 1 * zeros(3,3) * [1 1 0]' = [1 0 1]';

% learning rule
%W(t+1) = W(t) + ai(t+1) * Aj(t);
% example
%zeros(3,3) + [1 0 1]' * [1 1 0]



























