function [] = hasselmo9402_working_manip()
% hasselmo_2002 theta model

%% runtime parameters
p.nTrls = 1;
p.nTSteps = 25;
p.nCA1cells = 2;
p.nCA3cells = 3;
p.nECcells = 2;

p.dt = 0.005; % Default = 5 ms
p.thF = 8; % Default = 8 Hz
p.nCycles = 6; % How many cycles would we like per trial. Default = 6? 
p.stepsPerCycle = ceil(((1/p.thF)/p.dt));
p.phaseStep = (2*pi)/p.stepsPerCycle; % number of radians to increment theta phase by for each time steps
p.k = .5; % Max synaptic growth. Default = .5
p.lrate = 0.01; % learning rate, per Ehren's experiment
p.thetaScale = 1; % X from Hasselmo et al. 2002, p 799
% [0 to 1] smaller means theta has less effect on functional strength
% paper is vague about whether different values are used for the different
% layers (CA3 vs EC), for now we'll assume that it is a single fixed
% parameter. Default = 1?
p.time = linspace(0,1,p.nCycles*p.stepsPerCycle);


%% 
slopes = [.25, 1, 5];

CAones = nan(2,p.nCycles*p.stepsPerCycle,length(slopes));
achs = nan(1,p.nCycles*p.stepsPerCycle,length(slopes));
for i = 1:length(slopes)
    p.slope = slopes(i);
    [a1, pha] = switchSlope(p);
    CAones(:,:,i) = a1.CA1;
    achs(:,:,i) = pha.AChLvls;
    
end

figure; 
cnt = 1;
for ii = 1:2:length(slopes)*2
subplot(3,2,ii);
plot(p.time, CAones(:,:,cnt)', 'LineWidth',3);
ylabel('a_{CA1}\it(t) \rm(AU)','FontSize',14);
if cnt == 3, xlabel('Time (s)','FontSize',14); end
yticks([-.2 0 .5]); xticks([0 .5 1]);
ax = gca; ax.FontSize = 12;

subplot(3,2,ii+1);
plot(p.time, achs(:,:,cnt), 'LineWidth',3);
ylabel('\bf\psi\rm\it(t) \rm(AU)','FontSize',14);
text(.85,.25,['\gamma = ' num2str(slopes(cnt))], 'FontSize',12)
if cnt == 3, xlabel('Time (s)','FontSize',14); end
yticks([0 1]); xticks([0 .5 1]);
ax = gca; ax.FontSize = 12;

cnt = cnt + 1;
end



% subplot(3,2,3);
% plot(p.time, CAones(:,:,1)');
% 
% subplot(3,2,4);
% plot(p.time, achs(:,:,1));
% 
% subplot(3,2,5);
% plot(p.time, CAones(:,:,1)');
% 
% subplot(3,2,6);
% plot(p.time, achs(:,:,1));


end


function [a,pha] = switchSlope(p)

%todo: 
% add plots for ca3 and ec activity - done
% add ach to the model - done
% make a struct to hold data? - done

% can we get EC and CA3 activity separated?
% Learning? - simple vector forgeting? TMaze? 
% OLM switch? 
% show that hasselmo is wrong
% add the anon funct to other areas

dbstop if error;



% allocate arrays
a.CA1 = nan(p.nCA1cells,p.nTSteps);
a.CA3 = nan(p.nCA3cells,p.nTrls); 
a.EC  = nan(p.nECcells, p.nTrls);
w.CA3 = eye(p.nCA1cells, p.nCA3cells);
w.EC  = eye(p.nCA1cells, p.nECcells); % identity matrix p.801
tempXprod = nan(p.nCA1cells,p.nCA3cells,p.stepsPerCycle);


plts = 0;
%% Task 
stage = 1;

% for each trial, run the theta model 
for trl = 1:p.nTrls
  % initialize first timestep
  % NOTE: THESE NEED TO BE SET FOR EACH PHASE OF THE TASK, TRL
  a.CA3(:,1) = [1; 1; 0]; % rand(p.nCA3cells,1) > .5;
  a.EC(:,1)  = [1; 0]; % rand(p.nECcells, 1) > .5;
  
  % run the model
  [a, tempXprod, pha, theta, syn] = runTheta(a,w,tempXprod,p,stage);
  
  % compute weight updates after each theta cycle
  dw.CA3(:,:) = sum(tempXprod,3);
  
  fprintf('\nTrial %i\n',trl);
  if any(isnan(dw.CA3)), keyboard, end % Just in case something goes wrong...
  %w.CA3, dw.CA3
  w.CA3 = w.CA3 + p.lrate .* dw.CA3; %   w.CA3 = w.CA3 + dw.CA3;
  w.CA3 = min(w.CA3, p.k);
  %w.CA3
  
  %if trl == p.nTrls, plotStateVariables(theta,a,pha, trl); end
  
end




fprintf('Test initial learning\n');
a.CA3(:,1) = [0; 1; 0]; % rand(p.nCA3cells,1) > .5;
[a, tempXprod, pha, theta, syn] = runTheta(a,w,tempXprod,p,stage);
trl = 0;

if plts
    plotStateVariables(theta,a,pha, trl,p.time); 
    plotRegionAct(syn.EC3(1,:), syn.CA3(2,:),p.time);
end


if 0 % crosses may be indicative of the results found in colgin et al 2009
   figure; 
   subplot(2,1,1);
   hold on;
   plot(p.time,syn.EC3(1,:))
   plot(p.time,syn.CA3(2,:))
   
   subplot(2,1,2);
   plot(p.time,theta.LTP)
   
end

%keyboard;
end 


function [a,tempXprod, pha, theta, syn] = runTheta(a,w,tempXprod,p,stage)
  
  % allocate phase and theta params
  pha.EC = nan(p.nTSteps,1);   theta.EC = nan(p.nTSteps,1);
  pha.CA3 = nan(p.nTSteps,1);  theta.CA3 = nan(p.nTSteps,1);
  pha.LTP = nan(p.nTSteps,1);  theta.LTP = nan(p.nTSteps,1);
  
  pha.EC(1)  = 0;    theta.EC(1)  = (p.thetaScale/2) * sin(pha.EC(1))  + (1-(p.thetaScale/2));
  pha.CA3(1) = pi;   theta.CA3(1) = (p.thetaScale/2) * sin(pha.CA3(1)) + (1-(p.thetaScale/2));
  pha.LTP(1) = 0;    theta.LTP(1) = sin(pha.LTP(1));
  
  % initialize ACh activity
  AchPut = linspace(-10,10,p.nCycles*p.stepsPerCycle);
  pha.AChLvls = sigmoid(AchPut,p.slope);
  % .64 C_theta (equ 5) Cr = .8 Cl = ? 
  Cr = 1; Cl = .5; % These need to be updated to control for individual dynamics in CA3 vs EC 
  updateEC3 = @(ACh,Cl,theta,w,a) (ACh*Cl) * ((theta.*w) * a);
  updateCA3 = @(ACh,Cr,theta,w,a) (1 - ACh*Cr) * ((theta.*w) * a);
  
  syn.EC3(:,1) = updateEC3(pha.AChLvls(1),Cl,theta.EC(1),w.EC,a.EC(:,stage));
  syn.CA3(:,1) = updateCA3(pha.AChLvls(1),Cr,theta.CA3(1),w.CA3,a.CA3(:,stage));

  % initialize CA1 activity
  a.CA1(:,1) = syn.EC3(:,1) + syn.CA3(:,1); % eq 2.4 p.799

  tempXprod(:,:,1) = (theta.LTP(1) .* a.CA1(:,1)) * a.CA3(:,stage)';
  % run one theta cycle
  for t = 2:p.nCycles*p.stepsPerCycle
    
    % phasic input
    pha.EC(t)  = pha.EC(t-1)  + p.phaseStep;  theta.EC(t)  = (p.thetaScale/2) * sin(pha.EC(t))  + (1-(p.thetaScale/2));  % eq 2.2 p.799
    pha.CA3(t) = pha.CA3(t-1) + p.phaseStep;  theta.CA3(t) = (p.thetaScale/2) * sin(pha.CA3(t)) + (1-(p.thetaScale/2));  % eq 2.3 p.799
    pha.LTP(t) = pha.LTP(t-1) + p.phaseStep;  theta.LTP(t) = sin(pha.LTP(t));                                            % eq 2.5 p.799
    
    % after the new updating function was written
    syn.EC3(:,t) = (pha.AChLvls(t)*Cl) * ((theta.EC(t) .* w.EC) * a.EC(:,stage));
    syn.CA3(:,t) = (1 - pha.AChLvls(t)*Cr) * ((theta.CA3(t) .* w.CA3) * a.CA3(:,stage));
    
    a.CA1(:,t) = syn.EC3(:,t) + syn.CA3(:,t); % eq 2.4 p.799

    tempXprod(:,:,t) = (theta.LTP(t) .* a.CA1(:,t)) * a.CA3(:,stage)';
    
  end
  
  modEC = reshape(syn.EC3,2,[])';
  modCA3 = reshape(syn.CA3,2,[])';
  
  if 0
    keyboard;
    
    figure; hold on;
    time = linspace(0, 1000, size(reshape(syn.EC3,2,[])',1));
    plot(time,modEC);
    plot(time,modCA3);
    xlim([0,max(time)]);
    
    figure;
    plot(time, modEC+modCA3)
    
  end
  
end

function plotStateVariables_old(theta,a,pha,trl)
    figure;
    subplot(4,1,1);
    hold off; plot(theta.EC)
    hold on; plot(theta.CA3)
    ylabel('theta.CA3');
    legend('EC','CA3'); ylim([0 1]);
    %title(['End of training']);
    title(['trl ', num2str(trl)])
    
    subplot(4,1,2);
    plot(a.CA1');
    ylabel('act CA1');
    
    subplot(4,1,3);
    plot(theta.LTP)
    ylabel('theta LTP');
    
    subplot(4,1,4);
    plot(pha.AChLvls)
    ylabel('Ach');
end

function plotStateVariables(theta,a,pha,trl,t)
    figure;
    subplot(3,1,1); 
    plot(t,a.CA1','LineWidth',3);
    title('ACh-Theta Model Output','FontSize',16);
    ylabel('a_{CA1}\it(t) \rm(AU)','FontSize',14);
    yticks([-.2 0 .5]); xticks([0 .5 1]);
    ax = gca; ax.FontSize = 12;
    
    
    subplot(3,1,2);
    plot(t,pha.AChLvls,'LineWidth',3)
    ylabel('\bf\psi\rm\it(t) \rm(AU)','FontSize',14); 
    xlabel('Time (s)','FontSize',14);
    yticks([0 1]); xticks([0 .5 1]);
    ax = gca; ax.FontSize = 12;
    
    subplot(3,1,3);
    plot(t,theta.LTP,'LineWidth',3)
    ylabel('\theta_{LTP} \rm(AU)','FontSize',14);
    yticks([-1 0 1]); xticks([0 .5 1]);
    ax = gca; ax.FontSize = 12;
    
end

function plotRegionAct(EC, CA3, t)

figure; 
subplot(1,2,1);
plot(t,EC, 'r','LineWidth',3)
title('EC Activity', 'FontSize',16)
xlabel('Time (s)','FontSize',14); ylabel('Activity (AU)', 'FontSize',14);
yticks([0 .5]); xticks([0 .5 1]);
ax = gca; ax.FontSize = 12; 



subplot(1,2,2);
plot(t,CA3, 'b','LineWidth',3)
title('CA3 Activity', 'FontSize',16)
xlabel('Time (s)','FontSize',14); ylabel('Activity (AU)', 'FontSize',14);
yticks([0 .5]); xticks([0 .5 1]);
ax = gca; ax.FontSize = 12; 

end

function y = sigmoid(x,slope)
% slope > 1 steep, slope < 1 shallow
y = 1.0 ./ (1 + exp(slope*-x));
end
%%
