% hasselmo_2002 theta model

t = 5; 
cells = 5; 
a_CA1 = rand(cells, t, 1);
a_CA3 = rand(cells, t, 1) > .5; 
a_EC = rand(cells, t, 1) > .5;
w_CA3 = rand(cells, t, 1);
w_EC = rand(cells, t, 1);

% eq 2.1
% a_CA1 = w_EC * a_EC + w_CA3* a_CA3;

X = rand(cells, t, 1); % magnitude of synaptic currents of layers: 0 < X < 1 
phase_EC = pi/2;  
phase_CA3 = 3*pi/2;  % these need to be offset

% phasic input
theta_EC = X / 2 * sin(t + phase_EC) + (1 - X / 2); 
theta_CA3 = X / 2 * sin(t + phase_CA3) + (1 - X / 2);

% eq 2.4
a_CA1 = (theta_EC * w_EC * a_EC) + (theta_CA3 * w_CA3 * a_CA3);

plot(a_CA1)
