%% flineseek
% Linewise file seeking function.
%
% outpos = flineseek(fid,lineID)
%
% This function moves the in-file pointer to the beginning of the specified
% line. The file must already be open for reading. It returns the new position.
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master Linewise File Seeking Function code

function outpos = flineseek(fid,lineID)
% The purpose of this function is to move the in-file pointer to the beginning
% of the specified line, and return the new position.
%
% INPUTS/OUTPUTS
% fid - File ID to seek in.
% lineID - Line number of the line to move the pointer to the beginning of.
% outpos - New pointer position.

% Moving the pointer to the beginning
frewind(fid);

% Seeking through the file line-by-line
for lid = 1:(lineID-1)
   fgets(fid);
end;

% Retrieving the new pointer
outpos = ftell(fid);
end