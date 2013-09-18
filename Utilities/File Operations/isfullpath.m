%% isfullpath
% Relative/full path checking function.
%
% iffull = isfullpath(rfpath)
%
% This function checks if one or more file/directory path(s) is/are full
% (absolute) or relative. Its algorithm is listed in the Algorithm cell below.
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Algorithm
%
% On a Windows machine:
%  A path is considered full if it contains the ':' symbol as the second symbol
%  of the path string. In all other cases it is necessarily relative. This does
%  not necessarily imply that it is legal.
%
% On a Linux/Unix/Mac machine:
%  A path is considered full if it starts with the '/' symbol. In all other
%  cases it is necessarily relative. This does not necessarily imply that it is
%  legal.

%% Main path fullness checking function

function iffull = isfullpath(rfpath)
% The purpose of this function is to check if one or more file/directory path(s)
% is/are full(absolute) or relative.
%
% INPUTS/OUTPUTS
% rfpath - File/directory path(s) to check if it/they are full. String OR cell
%          array of strings.
% iffull - Whether or not rfpath is/are full. Boolean scalar OR array the size
%          of rfpath.

% Forcing cell array representation
if ischar(rfpath)
   rfpath = {rfpath};
elseif iscellstr(rfpath)
   % Already a cell array of strings
else
   error('Unexpected path format!!!');
end;

% Checking platform
if ispc
   isabs = @(rfpath) ischar(rfpath) && (length(rfpath) > 2) && ...
      (rfpath(2) == ':');
elseif isunix || ismac
   isabs = @(rfpath) ischar(rfpath) && (length(rfpath) >= 1) && ...
      (rfpath(1) == '/');
end;

% Checking path fullness
iffull = cellfun(isabs,rfpath);
end