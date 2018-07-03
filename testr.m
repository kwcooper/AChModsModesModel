function testr()



function [act] = g(A, Thresh)
  act = max(A - Thresh, 0);
end

% Test the g function -> outdated
num = -10:10;
sm = ones(size(num)) ./ num;
sm(find(isinf(sm))) = 0; 

for i = 1:size(sm,2)
  disp(sm(i))
  gO(:,i) = g(sm(i), .4);
end
figure; plot(sm,gO); title('G(CA1) | Thresh = 0.4')


xi = 2.0;
V = 3.0;
A_CA1 = [-1 -1 -1 -1]; % [.2 .4 0 .1];
%psi = (1 + exp(xi * sum(g(A_CA1,0.4) - V)))^-1;
psi = (1 + exp(xi * sum(max(A_CA1 - 0.4, 0) - V)))^-1;


% test the psi function 
num = -10:10;
for i = 1:size(num, 2)
  xi = 2.0;
  V = 3.0;
  A_CA1 = num(i); %ones(1,num(i)); % [.2 .4 0 .1];
  %psi = (1 + exp(xi * sum(g(A_CA1,0.4) - V)))^-1;
  psi = (1 + exp(xi * sum(max(A_CA1 - 0.4, 0) - V)))^-1;
  psiO(i) = psi;
end
figure; plot(num,psiO); title('Psi Output');

% now add the c's to the show
% Hasselmo altered these extensivly in a full parameter search
% can use figure 7 in the 1994 paper for optimal value ranges

% example from the paper for testing
C_theta = 0.64;
C_H = 0.8;
C_R = 0.8;

xi = 2.0;
V = 3.0;
A_CA1 = [0 0 0]; % [.2 .4 0 .1];
%psi = (1 + exp(xi * sum(g(A_CA1,0.4) - V)))^-1;
psi = (1 + exp(xi * sum(max(A_CA1 - 0.4, 0) - V)))^-1;

synTran = 1 - psi * C_R;


end

