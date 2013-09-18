%% repvectelems
% Vector element replication function.
%
% V = repvectelems(v,enums)
%
% This function takes a vector of objects and repeats each element a specified
% number of times. The repeated objects are all concatenated into a new vector
% with the same orientation as the input vector (default: horizontal) and
% returned.
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Vector element replication function

function V = repvectelems(v,enums)
% The purpose of this function is to replicate each element in a vector a given
% number of times.
%
% INPUTS/OUTPUTS
% v - Vector of elements to repeat.
% enums - Vector of numbers of repetitions for each element in v.
% V - Concatenated vector of repeated elements.

% Checking emptiness
if isempty(v)
   V = v;
   return;
end;

% v is nonempty
% Checking sizing
assert(length(v) == length(enums),'Mismatched vector lengths!!!');
N = length(v);

% Identifying orientation
if iscolumn(v) && ~isscalar(v)
   ifvert = true;
else
   ifvert = false;
end;

% Replicating (using repmat when possible, indexing otherwise)
if isscalar(v)
   V = repmat(v,1,enums);
else
   subs = cell(1,N);
   for vid = 1:N
      subs{vid} = vid * ones(1,enums(vid));
   end;
   V = v([subs{:}]);
end;

% Checking orientation
if ifvert
   V = V';
end;
end