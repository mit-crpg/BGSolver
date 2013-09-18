%% runAllmatFileScripts
% Script for running all .mat file building scripts in the mat_Files subpackage.
%
% This script runs the scripts listed in the Running .mat file scripts cell. It
% also provides output for the user.
%
% Package:    BGSolver v1.03
% Subpackage: mat_Files
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Running .mat file scripts

% BGSD_defaults.mat
BGSD_defaultsBuild;
fprintf('Created the following .mat file(s):\n');
fprintf('- BGSD_defaults.mat\n');

% ti_info.mat
ti_infoBuild;
fprintf('- ti_info.mat\n');