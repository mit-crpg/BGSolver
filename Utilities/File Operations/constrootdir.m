%% constrootdir
% Absolute root directory path construction function.
%
% rootdir = constrootdir(fname)
%
% This function constructs the absolute root directory path for a given file. If
% only the file name is given, without an absolute or a relative path, the
% current directory will be given as the absolute root directory path.
%
% Package:    BGSolver v1.03
% Subpackage: Utilities
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master Root Directory Path Construction Function code

function rootdir = constrootdir(fname)
% The purpose of this function is to construct the absolute root directory path
% for a given file.
%
% INPUTS/OUTPUTS
% fname - File name to construct the absolute root directory path for.
% rootdir - Absolute root directory path.

fname = GetFullPath(fname);
rootdir = fileparts(fname);
end