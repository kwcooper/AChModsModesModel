
% % sigmoid
% beta = .1; % changes steepness of the curve
% theta = .5;
% a = -100:100;
% %sigmoid = 1+exp(beta * (a - theta)^ -1);
% 
% out = [];
% for i = 1:size(a,2)
%   out(i) = 1 / (1+exp(-beta * (a(i) - theta)));
% end
% 
% figure; plot(out); title
% 
% 
% a = -100:100;
% s = 0;
% out = [];
% for i = 1:size(a,2)
%   out(i) = exp(-a(i) * s)/s;
% end
% 
% figure; plot(out); title


% phase precession model
t = linspace(0,2.5,200);
dend = [];
for i = 1:size(t,2)
  dend(i) = cos(2*pi*6*t(i));
end

som = [];
for i = 1:size(t,2)
  som(i) = cos(2*pi*6.4*t(i));
end

sum = [];
for i = 1:size(t,2)
  sum(i) = cos(2*pi*6*t(i)) + cos(2*pi*6.4*t(i));
end

thresh = 1.5;
spikes = sum > thresh;

figure;
subplot(4,1,1)
plot(dend); title('Dendrite');
subplot(4,1,2)
plot(som); title('Soma');
subplot(4,1,3)
plot(sum); title('V(t): Sum of Oscillations');
subplot(4,1,4)
plot(spikes); title('Spikes');
