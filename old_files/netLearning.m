% learning in networks


% learning rule example
s1 = [0 1 1 0];
s2 = [0 1 0 1];

WM = s1' * s2; WM

% activation rule
s2p = [0 0 0 1];
a1 = WM * s2p'; a1

%new learning; new network
s3 = [1 0 1 0];
s4 = [1 0 0 1];

WM2 = s3' * s4; WM2

% activation rule
s4p = [0 0 1 0];
a2 = WM2 * s4p'; a2
a22 = WM2 * s4p'; a22


% test on larger matricies
n = 10;

% create input
p1 = randi([0 1], n,1);
p2 = randi([0 1], n,1);

p1 * p2'



