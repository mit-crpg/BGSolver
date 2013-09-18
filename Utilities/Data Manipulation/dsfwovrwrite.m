%% dsfwovrwrite
% Data structure fieldwise overwrite function.
%
% dsbase = dsfwovrwrite(dsbase,dsovr)
%
% This function takes a data structure and overwrites its fields with fields
% from another data structure, while not removing the fields the original data
% structure had that the new data structure does not.
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Data structure fieldwise overwrite function code

function dsbase = dsfwovrwrite(dsbase,dsovr)
% The purpose of this function is to overwrite fields in a base data structure
% using fields with the same names from another data structure, or add new
% fields if the base data structure doesn't have them.
%
% INPUTS/OUTPUTS
% dsbase - Base data structure.
% dsovr - Overwriting data structure.

fldovr = fieldnames(dsovr);

for fldid = 1:length(fldovr)
   dsbase.(fldovr{fldid}) = dsovr.(fldovr{fldid});
end;
end