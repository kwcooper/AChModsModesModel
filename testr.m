function testr()



function [act] = g(A, Thresh)

if A - Thresh > 0
  act = A - Thresh;
elseif A - Thresh >= 0
  act = 0;
else
  act = 0;
end

end

% Test the g function
num = -10:10;
sm = ones(size(num)) ./ num;
sm(find(isinf(sm))) = 0; 

for i = 1:size(sm,2)
  disp(sm(i))
  gO(:,i) = g(sm(i), .4);
end
figure; plot(gO); title('G(CA1) | Thresh = 0.4')


xi = 2.0;
%psi = 1 + exp(xi);


end