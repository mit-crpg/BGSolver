%% rmfieldcheck
% Pre-checking field removal function.
%
% ds = rmfieldcheck(ds,dsfldnames)
%  OR
% ds = rmfieldcheck(ds,dsfldname)
%
% This function checks whether or not a data structure has a field(s) with a
% given name(s), and if it does - removes it (them).
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Pre-checking field removal function

function ds = rmfieldcheck(ds,dsfldnames)
% The purpose of this function is to check whether or not a data structure has a
% field(s) with a given name(s), and if it does - to remove it (them).
%
% INPUTS/OUTPUTS
% ds - Data structure.
% dsfldnames - Field name(s) to check and remove, if present. Either cell array
%              of strings, or a string.

if iscellstr(dsfldnames)
   ids = isfield(ds,dsfldnames);
   ds = rmfield(ds,dsfldnames(ids));
elseif ischar(dsfldnames)
   if isfield(ds,dsfldnames)
      ds = rmfield(ds,dsfldnames);
   end;
else
   error('dsfldnames is of invalid data type!!!');
end;
end