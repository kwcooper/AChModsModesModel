

% task

K = 1;

%% Learn left arm (n)
a.EC.n  = [1 0];
a.CA3.n = [1 .5 0]; %[1 1 0];
w.EC.n = eye(2, 2);
w.CA3.n = K * a.CA3.n'* a.EC.n; % eq (2.8) % 0 .* ones(3,2);

% Activation % CA1 = [1 0]
a.CA1.n = a.EC.n * w.EC.n + a.CA3.n * w.CA3.n; % eq (2.1)

% does a.CA1.n == a.EC.n? 
 
%% Error (e)
a.EC.e  = [0 0];
a.CA3.e = [.5 1 0];
w.EC.e = eye(2, 2);
w.CA3.e = w.CA3.n +  a.CA3.n' .* a.CA1.n;    

% Activation % CA1 = [0 1]
a.CA1.e = a.EC.e * w.EC.e + a.CA3.e * w.CA3.e;


%% Correct Reversal (c)
a.EC.c  = [0 1];
a.CA3.c = [0 1 1];
w.EC.c = eye(2, 2);
w.CA3.c = w.CA3.e +  a.CA3.e' .* a.CA1.e;  

% Activation % CA1 = [0 1]
a.CA1.c = a.EC.c * w.EC.c + a.CA3.c * w.CA3.c;


%% Test (r)
a.EC.r  = [0 0];
a.CA3.r = [1 1 1];
w.EC.r = eye(2, 2);
w.CA3.r = w.CA3.c +  a.CA3.c' .* a.CA1.c;

% Activation % CA1 = [0 1]
a.CA1.r = a.EC.r * w.EC.r + a.CA3.r * w.CA3.r;


%% Scratchpad

if 0
% pg 802 2. 
a.CA3.n * a.CA3.e';     % should equal 1
[1 .5 0] * [.5 1 0]';   %

% pg 802 3. 
a.CA3.n * a.CA3.c';     % should equal 0 
[.5 .5 0] * [0 .5 .5]';
[0 .25 .25] * [.25 .25 0]';
a.CA3.e * a.CA3.c';     % should equal 0 
a.EC.c * a.EC.n';       % should equal 0
a.EC.c * a.EC.c';       % should equal 1
a.EC.n * a.EC.n';       % should equal 1


% pg 804 2.4
a.CA3.n * a.CA3.r;      % should equal 1
a.CA3.c * a.CA3.r;      % should equal 1
%[q1 q2 0] * [.5 1 0]' == 1 && [] * []' == 0

end






























