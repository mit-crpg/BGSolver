%% Sample_RLC_Build
% .mat file building script for the Series RLC circuit BGSD in the
% Sample_BGSDs subpackage.
%
% This script builds the .mat file which contains runopts, bgs and options
% data structures for a fully linear series RLC circuit with a sinusoidal
% voltage source.
%
% Package:    BGSolver v1.03
% Subpackage: Sample_BGSDs
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Setting parameters

C = 0.5; % Capacitance 0.5 F
L = 1.5; % Inductance 1.5 H
R = 100; % Resistance 100 Ohm
Va = 10; % Voltage source amplitude 10 V
f = 10; % Voltage source frequency 10 Hz

omega = 2*pi/f;
V = @(t) Va * sin(omega*t); % Voltage source function

Vc = @(q) (1/C)*q; % Capacitor-enforced voltage function

fnameData = 'RLC_Data.mat';
fnameBGSD = 'RLC_BGSD.mat';

%% Storing supporting quantities in a .mat file

save(fnameData,'Va','omega','C');

%% Constructing runopts

% Runtime and Recording Options
runopts.verbose = 2;
% Time Integrator Options
runopts.iispec = 0;
runopts.ii = [0,100];
runopts.addtio = 3;
runopts.addtiolist = {'MaxOrder','2';
                      'MaxStep','0.01';
                      'NonNegative','[]'};
runopts.tidisp = 1;
% Linked Files and Directories
runopts.lflist = {fnameData};
% Post-Processor Options
runopts.epp = 2;
runopts.eppIDs = [4,2];
runopts.ppdisp = 1;

%% Constructing bgs

% Bond Connectivity Map
bgs.bonds = uint16([2,1;
                    3,2;
                    4,2;
                    5,2]);

% Element List
bgs.elements = preallocEL(5);
                        
bgs.elements(1).etype = 'MSE';
bgs.elements(1).exprtype = 'NMC';
bgs.elements(1).mvars = sym('t');
bgs.elements(1).expr = V;

bgs.elements(2).etype = '1';

bgs.elements(3).etype = 'C';
bgs.elements(3).exprtype = 'NE';
bgs.elements(3).x0 = 0;
bgs.elements(3).expr = Vc;
bgs.elements(3).phx = 1;
bgs.elements(3).nn = true;

bgs.elements(4).etype = 'I';
bgs.elements(4).exprtype = 'CC';
bgs.elements(4).x0 = 0;
bgs.elements(4).expr = 1/L;
bgs.elements(4).phx = 2;
bgs.elements(4).nn = true;

bgs.elements(5).etype = 'R';
bgs.elements(5).exprtype = 'CC';
bgs.elements(5).expr = R;

%% Storing runopts and bgs

save(fnameBGSD,'runopts','bgs');

%% Clearing memory

clear('C','L','R','Va','f','omega','V','Vc','fnameData','fnameBGSD',...
   'runopts','bgs');