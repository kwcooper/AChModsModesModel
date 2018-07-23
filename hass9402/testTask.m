% test task 


%       CA3      EC
% T1:  1 1 0    1 0
% T2:  0 1 1    0 1
% Ch:  0 0 1     ?    should be 0 1  not 1 0... 
% Ch2: 0 1 0     ?    Will it be 1 1?



a.CA3(:,1) = [0;1;1]; 
a.EC(:,1)  = [0;1]; 

%a.CA3(:,2) = [0;1;1]; 
%a.EC(:,2)  = [0;1]; 

%a.CA3(:,3) = [0;0;1;]; 


% test activation 
w = EC(:,1) * CA3(:,1);


