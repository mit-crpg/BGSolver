%% getvtopts
% Variable type options retrieving function.
%
% vtopts = getvtopts(v,vars2check,[vars2order])
%
% This function constructs a variable type options data structure for getvvids
% and getevids functions. The variable type options data structure contains
% information about whether or not to check the variable, whether or not the
% variable is ordered, etc.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       January 15, 2013
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master variable type options retrieving function code

function vtopts = getvtopts(v,vars2check,vars2order)
% The purpose of this function is to retrieve the options (whether or not to
% check, whether or not the variable type is ordered, etc.) of a given variable
% type based on inputs.
%
% INPUTS/OUTPUTS
% v - Variable type character. REQUIRED.
% vars2check - Variable types to extract. String. REQUIRED.
% vars2order - Order tracking status vector. Boolean vector. OPTIONAL. By
%              default, order is not tracked.
% vtopts - Variable type options data structure. Its fields are:
%  .ifcheck - Whether or not to check the corresponding variable type. If this
%             field is false, the other fields are empty.
%  .idInOut - Variable type's index in output cell. uint8.
%  .iforder - Whether or not the corresponding variable type is ordered.
%             Boolean.

[vtopts.ifcheck,vtopts.idInOut] = ismember(v,vars2check);
vtopts.idInOut = uint8(vtopts.idInOut);

vtypnum = uint8(length(vars2check));

vtopts.iforder = false;
if nargin == 3 && vtopts.ifcheck && vtopts.idInOut <= vtypnum && ...
      vtypnum <= length(vars2order)
   vtopts.iforder = vars2order(vtopts.idInOut);
end;
if v == 't'
   vtopts.iforder = false;
end;
end