% hasselmo_2002 theta model

duration = 10; 
cells = 5; 

a_CA1 = nan(cells, duration, 1);
a_CA3 = nan(cells, duration, 1) > .5; 
a_EC = nan(cells, duration, 1) > .5;
w_CA3 = nan(cells, duration, 1);
w_EC = nan(cells, duration, 1);

a_CA1(:,1) = rand(cells, 1, 1);
a_CA3(:,1) = rand(cells, 1, 1) > .5; 
a_EC(:,1) = rand(cells, 1, 1) > .5;
w_CA3(:,1) = rand(cells, 1, 1);
w_EC(:,1) = rand(cells, 1, 1);

% eq 2.1
%a_CA1 = w_EC .* a_EC + w_CA3 .* a_CA3;

%X = rand(cells, t, 1); % magnitude of synaptic currents of layers: 0 < X < 1 

phase_EC = pi/2;  
phase_CA3 = 3*pi/2;  % these need to be offset

for t = 2:duration 
X = rand(cells, 1, 1);
% phasic input
theta_EC = X / 2 * sin(t + phase_EC) + (1 - X / 2); 
theta_CA3 = X / 2 * sin(t + phase_CA3) + (1 - X / 2);

% eq 2.4
a_CA1(:,t) = (theta_EC .* w_EC(:,t-1) .* a_EC(:,t-1)) + (theta_CA3 .* w_CA3(:,t-1) .* a_CA3(:,t-1));

w_CA3(:,t) = rand(cells, 1, 1);
w_EC(:,t) = rand(cells, 1, 1);
end

plot(a_CA1')

