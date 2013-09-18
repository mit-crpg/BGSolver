%% cell2str
% Cell vector of strings to string conversion function.
%
% str = cell2str(cellofstrs,[delim])
%
% This function converts a cell vector of strings into a single string,
% separated by an optional delimiter. If a delimiter is not specified, a
% whitespace is used. The last entry is not followed by a delimiter.
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Cell Vector of Strings to String Conversion Function Code

function str = cell2str(cellofstrs,delim)
% The purpose of this function is to convert a cell vector of strings into a
% single string.
%
% INPUTS/OUTPUTS
% cellofstrs - Cell vector of strings to convert.
% delim - Delimiter to separate strings with in the resulting string. OPTIONAL.
%         Default value: ' '.
% str - Resulting string.

% Identifying delimiter
if nargin == 1
   delim = ' ';
end;

% Checking if there is more than one string
if isscalar(cellofstrs)
   str = cellofstrs{1};
else
   % Constructing format
   fmt = ['%s' delim];
   
   % Converting the string
   str = [sprintf(fmt,cellofstrs{1:(end-1)}), cellofstrs{end}];
end;
end