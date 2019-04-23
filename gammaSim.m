% Gamma LFP Model

function gammaSim() 

fs = 600;              % Sampling frequency (samples per second)
dt = 1/fs;             % seconds per sample
StopTime = 1;        % seconds
t = (0:dt:StopTime)';  % seconds


F = 8;                 % theta frequency (hertz)
theta = sin(2*pi*F*t);

F = 55; % Hz                 
sGamma = sin(2*pi*F*t);
sgTheta = -sin(2*pi*8*t+pi/2);

F = 80; % Hz                
fGamma = sin(2*pi*F*t);

gScale = .3;
signal = theta + gScale * (fGamma .* relu(theta)) + gScale * (sGamma .* relu(sgTheta));

%%
figure;
numPlts = 4;
subplot(numPlts,1,1)
plot(t,theta); hold on; plot(t,sgTheta)
ylabel('\theta''s','fontsize',14)

subplot(numPlts,1,2)
plot(t,gScale * (sGamma .* relu(-sin(2*pi*8*t+pi/2))))
ylabel('Slow \gamma (CA3)','fontsize',11)

subplot(numPlts,1,3)
plot(t,gScale * (fGamma .* relu(theta)))
ylabel('Fast \gamma (EC)','fontsize',11)

subplot(numPlts,1,4)
%plot(theta + gScale * (fGamma .* relu(theta)) + gScale * (sGamma .* relu(-theta)))
plot(t,signal)
ylabel('CA1 LFP','fontsize',11)
xlabel('Time (S)','fontsize',11)


%% Generate theta mod gamma figure

[modindex, thetarange, gammarange, powPhsDists, bincenters,thetaamps_M,gammaamps_M,stdVals] = thetaModGamma_nonCMB(signal, fs,'filtParams',2,'stdGamma',1,'gammarange',[10:2:140],'thetarange',8);
%thetaPhaseAxis = [linspace(90,360,18)];
thetaPhaseAxis = linspace(-180,180,36); 
figure; imagesc(thetaPhaseAxis,gammarange,powPhsDists)
colormap 'jet';
title('Simulated theta mod gamma')
ylabel('Frequency (Hz)'); xlabel('\theta Phase')
%xticks([-80 90 270]);
%xticks(thetaPhaseAxis)
%
%h = colorbar; set(get(h,'label'),'string','Power (\mew V^{2})'); 
h = colorbar; ylabel(h, 'Power (\muV^{2})','rotation', 270, 'position', [2.7, 1.6, 0], 'fontsize',11);
axis xy; 

%% Square wave stuff
if 0
  sqwave = square(4*pi*t,25);
  figure; plot(t,sqwave);
end


end

function y = relu(x)
y = max(0,x);
end

function y = sig(x)
y = 1./(1 + exp(-x));
end