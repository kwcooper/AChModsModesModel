% test task 


%       CA3      EC
% T1:  1 1 0    1 0
% T2:  0 1 1    0 1
% Ch:  0 0 1     ?    should be 0 1  not 1 0... yes
% Ch2: 0 1 0     ?    Will it be 1 1? yes


%%

w(:,:,1) = zeros(2,3);

% Trial 1
CA3(:,1) = [0;1;1]; 
EC(:,1)  = [0;1];  

% learning  
w(:,:,2) = w(:,:,1) + (EC(:,1) * CA3(:,1)');

% Trial 1
CA3(:,2) = [1;1;0]; 
EC(:,2)  = [1;0]; 

% learning  
w(:,:,3) = w(:,:,2) + (EC(:,2) * CA3(:,2)');

% Checking
CA3(:,2) = [1;0;0];
EC(:,2) = [0;0];
CA3(:,3) = [0;0;1];
EC(:,3) = [0;0];
CA3(:,4) = [0;1;0];
EC(:,4) = [0;0];
% recall activation
EC(:,2) + w(:,:,3) * CA3(:,2);
EC(:,3) + w(:,:,3) * CA3(:,3);
EC(:,4) + w(:,:,3) * CA3(:,4);

%%
% Test accuracy measurement
a = [3 1 1];
b = [0 0 0];
%sqrt(1^2 + 1^2 + 0^2) hardcode norm

cosSim(a,b)
