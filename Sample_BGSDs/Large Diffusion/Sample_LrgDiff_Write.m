%% Sample_LrgDiff_Write
% BGSD file writing script for the large diffusion BGSD in the
% Sample_BGSDs subpackage.
%
% This script reads the .mat file which contains the runopts and bgs data
% structures, completes them, and writes the outcome in a BGSD file.
%
% Package:    BGSolver v1.03
% Subpackage: Sample_BGSDs
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Loading data structures

load('LrgDiff_BGSD.mat');

%% Constructing options

opmode = 'writeBGSD';

%% Completing the BGSD data structure

[runopts,bgs] = completeBGSDdata(runopts,bgs,'opmode',opmode);

%% Writing the BGSD file

writeBGSD(runopts,bgs);