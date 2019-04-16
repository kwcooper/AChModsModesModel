

% task

% Learn left arm
a.EC_1  = [1 0];
a.CA3_1 = [1 1 0];
w.EC_1 = eye(2, 2);
w.CA3_1 = 0 .* ones(3,2);

% Activation % CA1 = [1 0]
a.CA1_1 = a.EC_1 * w.EC_1 + a.CA3_1 * w.CA3_1;


% Error
a.EC_2  = [0 0];
a.CA3_2 = [1 1 0];
w.EC_2 = eye(2, 2);
w.CA3_2 = w.CA3_1 +  a.CA3_1' .* a.CA1_1;    

% Activation % CA1 = [0 1]
a.CA1_2 = a.EC_2 * w.EC_2 + a.CA3_2 * w.CA3_2;


% Correct Reversal
a.EC_3  = [0 1];
a.CA3_3 = [0 1 1];
w.EC_3 = eye(2, 2);
w.CA3_3 = w.CA3_2 +  a.CA3_2' .* a.CA1_2;  

% Activation % CA1 = [0 1]
a.CA1_3 = a.EC_3 * w.EC_3 + a.CA3_3 * w.CA3_3;


% Test
a.EC_4  = [0 0];
a.CA3_4 = [1 1 1];
w.EC_4 = eye(2, 2);
w.CA3_4 = w.CA3_3 +  a.CA3_3' .* a.CA1_3;

% Activation % CA1 = [0 1]
a.CA1_4 = a.EC_4 * w.EC_4 + a.CA3_4 * w.CA3_4;



