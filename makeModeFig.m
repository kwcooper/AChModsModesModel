% make gamma figures for thesis
function makeModeFig()

fs = 600;              % Sampling frequency (samples per second)
dt = 1/fs;             % seconds per sample
StopTime = 1;          % seconds
t = (0:dt:StopTime)';  % seconds

F = 8;                 % theta frequency (hertz)
theta = sin(2*pi*F*t);

F = 1;                 % ach frequency (need to look at data for this)
ach = sin(2*pi*F*t);

% ach & theta
if 0 figure; hold on; plot(t,theta); plot(t,ach); end
% ach mod theta
if 0 figure; plot(abs(theta + ach)); end


%%
% For one cycle get time period
T = 1/F ;
% time step for one time period
tt = 0:dt:T+dt;
d = sin(2*pi*F*tt);
%figure; plot(tt,d);


%% Gammas? 

F = 8;                 % theta frequency (hertz)
theta = sin(2*pi*F*t);

F = 60; % Hz                 
sGamma = sin(2*pi*F*t);
sgTheta = -sin(2*pi*8*t+pi/2);

F = 80; % Hz                
fGamma = sin(2*pi*F*t);

gScale = .3;

figure;
numPlts = 4;
subplot(numPlts,1,1)
plot(theta); hold on; plot(sgTheta)
ylabel('\theta''s','fontsize',14)

subplot(numPlts,1,2)
plot(gScale * (sGamma .* relu(-sin(2*pi*8*t+pi/2))))
ylabel('Slow \gamma (CA3)','fontsize',11)

subplot(numPlts,1,3)
plot(gScale * (fGamma .* relu(theta)))
ylabel('Fast \gamma (EC)','fontsize',11)

subplot(numPlts,1,4)
%plot(theta + gScale * (fGamma .* relu(theta)) + gScale * (sGamma .* relu(-theta)))
plot(theta + gScale * (fGamma .* relu(theta)) + gScale * (sGamma .* relu(sgTheta)))
ylabel('CA1 LFP', 'fontsize',11)

end

function y = relu(x)
y = max(0,x);
end

function y = sig(x)
y = 1./(1 + exp(-x));
end







