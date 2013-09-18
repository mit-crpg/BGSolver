%% srBGSolve
% Bond graph processing restart script.
%
% This script requests a Sorted .mat Results (SMR) file from the user using a
% dialog window, and restarts its bond graph processing using the frBGSolve. The
% outcome is stored in the T, X, E and F variables.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Retrieving SMR file from the user

mfname = openfiledlg({'*.mat','Sorted .mat Results (*.mat)'},...
   'Open an SMR file');

%% Processing the SMR file

[T,X,E,F] = frBGSolve(mfname);