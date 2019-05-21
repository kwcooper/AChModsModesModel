% Gamma LFP Model

function gammaSim3() 

fs = 600;              % Sampling frequency (samples per second)
dt = 1/fs;             % seconds per sample
StopTime = 1;        % seconds (~100s to run successfully)
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
figure('Position', [1 1 945 375]);
numPlts = 3;
subplot(numPlts,1,1)
plot(t,gScale * (sGamma .* relu(-sin(2*pi*8*t+pi/2))))
ylabel('Slow \gamma (CA3)','fontsize',11)
xticks([0 1]);
title('Simulated LFP (intra-cycle)')

subplot(numPlts,1,2)
plot(t,gScale * (fGamma .* relu(theta)))
ylabel('Fast \gamma (EC)','fontsize',11)
xticks([0 1]);


subplot(numPlts,1,3)
%plot(theta + gScale * (fGamma .* relu(theta)) + gScale * (sGamma .* relu(-theta)))
plot(t,signal)
ylabel('CA1 LFP','fontsize',11)
xlabel('Time (S)','fontsize',11)
xticks([0 1]);

if numPlts == 4
  subplot(numPlts,1,4)
  plot(t,theta); hold on; plot(t,sgTheta)
  ylabel('\theta''s','fontsize',14)
end


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


%% Thresholded gamma (Modes)



rt = relu(theta);
fgma = (fGamma .* rt);
zrt = rt ~= 0;
bw = buckWave(zrt);
%figure; plot(fgma.*bw);
fgmabw = fgma.*bw;

rt2 = relu(-sin(2*pi*8*t+pi/2));
sgma = (sGamma .* rt2);
zrt2 = rt2 ~= 0;
bw2 = buckWave(zrt2);
%figure; plot(sgma.*bw2);
sgmabw = sgma.*(1 - bw2); % flip bw2 so that they are offset

gScale = .3;
signal2 = theta + gScale * (fgmabw .* relu(theta)) + gScale * (sgmabw .* relu(sgTheta));



figure('Position', [1 1 945 375]); 
%hold on;
numPlts = 3;
subplot(numPlts,1,1);
plot(t,(gScale * sgmabw));
ylabel('Slow \gamma (CA3)','fontsize',11);
xticks([0 1]);
title('Simulated LFP (modes)')

subplot(numPlts,1,2);
plot(t,(gScale * fgmabw));
ylabel('Fast \gamma (EC)','fontsize',11);
xticks([0 1]);

subplot(numPlts,1,3);
%plot(theta + gScale * (fGamma .* relu(theta)) + gScale * (sGamma .* relu(-theta)))
plot(t,signal2)
ylabel('CA1 LFP','fontsize',11)
xlabel('Time (S)','fontsize',11)
xticks([0 1]);

if numPlts == 4
  subplot(numPlts,1,4)
  plot(t,theta); hold on; plot(t,sgTheta);
  ylabel('\theta''s','fontsize',14);
end


%% Generate theta mod gamma figure

[modindex, thetarange, gammarange2, powPhsDists2, bincenters,thetaamps_M,gammaamps_M,stdVals] = thetaModGamma_nonCMB(signal2, fs,'filtParams',2,'stdGamma',1,'gammarange',[10:2:140],'thetarange',8);
%thetaPhaseAxis = [linspace(90,360,18)];
thetaPhaseAxis = linspace(-180,180,36); 
figure; imagesc(thetaPhaseAxis,gammarange2,powPhsDists2)
colormap 'jet';
title('Simulated theta mod gamma (modes)')
ylabel('Frequency (Hz)'); xlabel('\theta Phase')
%xticks([-80 90 270]);
%xticks(thetaPhaseAxis)
%
%h = colorbar; set(get(h,'label'),'string','Power (\mew V^{2})'); 
h = colorbar; ylabel(h, 'Power (\muV^{2})','rotation', 270, 'position', [2.7, 1.6, 0], 'fontsize',11);
axis xy; 

%% Cycle stuff

  
  

end

function y = relu(x)
y = max(0,x);
end

function buckedWave = buckWave(wave)
% Takes a positive square wave, bucks

zrt = wave;
inds1 = find(diff(zrt)==1);
inds2 = find(diff(zrt)==-1);

ind2 =  inds2(2:4:end,:);
ind1 =  inds1(1:4:end,:);

buckedWave = zrt;
for ind = 1:size(ind2,1)
  buckedWave(ind1(ind):ind2(ind))=0;
end

end

function sm = findStreaks(x,i)
vd = diff([0 x' 0]); %the 0s ensure that there's always a pair of [+1 -1] in the diff
starts = find(vd == i);
ends = find(vd == -i); 
[sm, idx] = max(ends-starts);
end


function breakWave(x,t)
%https://www.mathworks.com/matlabcentral/answers/168261-i-need-a-code-to-separate-each-cycle-of-sine-wave
mcs = x .* circshift(x, [0 -1]);
zxix = find(mcs <= 0);
zx = nan(size(1:2:size(zxix,2)-1));
for k1 = 1:2:size(zxix,2)-1
  zx(k1) = interp1(x(zxix(k1):zxix(k1)+1), t(zxix(k1):zxix(k1)+1), 0);
end

if 1
  figure(1)
  plot(t, x)
  hold on
  plot(zx, zeros(size(zx)), '+r', 'MarkerSize',10)
  hold off
  grid
end
end

function y = sig(x)
y = 1./(1 + exp(-x));
end