function [result] = cosSim(a,b)
% returns the cosine similarity between two vectors
% checks for zero vectors, returning 0 if they occur

if sum(a) == 0 || sum(b) == 0
  result = 0;
else
  result = (a * b') / (norm(a) * norm(b));
end

end

