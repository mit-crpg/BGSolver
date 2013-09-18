%% Sample_RLC_Run
% Series RLC circuit BGSD processing script in the Sample_BGSDs subpackage.
%
% This script runs the processor on the series RLC circuit problem with modified
% options.
%
% Package:    BGSolver v1.03
% Subpackage: Sample_BGSDs
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Setting parameters

ovrrunopts.tit = 'bdf3';
ovrrunopts.tst = 1;
ovrrunopts.ftsspect = 0;
ovrrunopts.fts = 200;
ovrrunopts.sdisp = 1;
ovrrunopts.jsdisp = 1;
fname = 'test.bgsd';

%% Running the problem

[T,X,E,F] = fBGSolve(fname,'ovrrunopts',ovrrunopts);