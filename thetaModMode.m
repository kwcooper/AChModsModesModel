% theta model



% initial learning
AEC = [1 0]; %postsynaptic input
ACA3 = [1 1 0]; %presynaptic input

W = AEC' * ACA3;

AEC = [0 1];
ACA3 = [0 1 1];


% init params
t = ; 

% Encoding   (pi/2):  PLTP = 1;  PEC = 1; PCA3 = 0;
% Retrieval (3*pi/2): PLTP = -1; PEC = 0; PCA3 = 1;
PLTP = 1; % oscilates between 1 @ pi/2 & -1 @ 3*pi/2
PEC  = 1; % oscilates between 1 @ pi/2 &  0 @ 3*pi/2
PCA3 = 0; % oscilates between 1 @ pi/2 &  0 @ 3*pi/2



% define theta functions
thEC = .5 * sin(t + PEC) + .5;
thCA3 = .5 * sin(t + PCA3) + .5;
thLTP = sin(t + PLTP);

% activation rule
ai = thEC * AEC' + thCA3 * W * ACA3';

% learning

W_n = W +  









