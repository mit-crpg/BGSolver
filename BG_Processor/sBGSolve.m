%% sBGSolve
% Bond graph processing script.
%
% This script requests a BGSD file from the user using a dialog window, and
% processes it using the fBGSolve function. The outcome is stored in the T,
% X, E and F variables.
%
% Package:    BGSolver v1.02
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Retrieving BGSD file from the user

fname = openfiledlg({'*.bgsd',...
   'Bond Graph System Descriptor (*.bgsd)'},'Open a BGSD file');

%% Processing the BGSD file

[T,X,E,F] = fBGSolve(fname);