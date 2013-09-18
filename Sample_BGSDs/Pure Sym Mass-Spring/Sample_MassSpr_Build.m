%% Sample_MassSpr_Build
% .mat file building script for the purely symbolic mass-spring BGSD in the
% Sample_BGSDs subpackage.
%
% This script builds the .mat file which contains runopts and bgs data
% structures for a purely symbolic (no numeric variables) mass-spring system
% with a constant source of effort (gravity).
%
% Package:    BGSolver v1.03
% Subpackage: Sample_BGSDs
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Setting parameters

k = 100; % Stiffness (N/m)
b = 0.2; % Dashpot friction coefficient (kg/s)
m = 1; % Mass (kg)
g = 9.8; % Gravitational acceleration (m/s^2)

fnameData = 'MassSpr_Data.mat';
fnameBGSD = 'MassSpr_BGSD.mat';

%% Storing supporting quantities in a .mat file

save(fnameData,'k','b','m','g');

%% Constructing runopts

% Runtime and Recording Options
runopts.verbose = 2;
runopts.log = 2;
% Time integrator options
runopts.ii = [0,10];
% Linked Files and Directories
runopts.lflist = {fnameData};

%% Constructing bgs

% Bond Connectivity Map
bgs.bonds = uint16([4,2;
                    3,2;
                    1,2;
                    2,5]);

% Element List
bgs.elements = preallocEL(5);
                        
bgs.elements(1).etype = 'C';
bgs.elements(1).exprtype = 'CC';
bgs.elements(1).x0 = 1;
bgs.elements(1).expr = k;

bgs.elements(2).etype = '1';

bgs.elements(3).etype = 'R';
bgs.elements(3).exprtype = 'CC';
bgs.elements(3).expr = b;

bgs.elements(4).etype = 'I';
bgs.elements(4).exprtype = 'CC';
bgs.elements(4).x0 = 0;
bgs.elements(4).expr = m;

bgs.elements(5).etype = 'SE';
bgs.elements(5).exprtype = 'CC';
bgs.elements(5).expr = m*g;

%% Storing runopts and bgs

save(fnameBGSD,'runopts','bgs');

%% Clearing memory

clear('k','b','m','g','fnameData','fnameBGSD','runopts','bgs');