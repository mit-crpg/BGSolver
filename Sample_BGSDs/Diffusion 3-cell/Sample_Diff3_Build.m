%% Sample_Diff3_Build
% .mat file building script for the 3-cell diffusion BGSD in the
% Sample_BGSDs subpackage.
%
% This script builds the .mat file which contains runopts, bgs and options
% data structures for a nonlinear 3-cell diffusion BGSD with linear storage
% but nonlinear, modulated and asymmetric diffusion elements.
%
% Package:    BGSolver v1.03
% Subpackage: Sample_BGSDs
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Setting parameters

C = 10; % A cell's thermal stiffness
Qa = 10; % Cell thermal sources' amplitudes
omega = 2*pi/0.5; % Angular frequency of heat source oscillations
R1 = 2; % Resistance of the left R2
R2 = 2.5; % Resistive amplitude of the right R2 at T2 = 0
R2gr = 0.1; % Rate of growth of the right R2's resistance with T2

Qf = @(t) Qa*sin(omega*t);
Q1 = @(Ts) (1/R1)*(Ts(1,:)-Ts(2,:)).^3;
Q2cond = @(T) 1./(R2+R2gr.*T);
Q2condMin = @(T) -1./(R2+R2gr.*T);

fnameData = 'Diff3_Data.mat';
fnameBGSD = 'Diff3_BGSD.mat';

%% Storing supporting quantities in a .mat file

save(fnameData,'Qa','omega','R1','R2','R2gr');

%% Constructing runopts

% Runtime and Recording Options
runopts.verbose = 1;
% Time Integrator Options
runopts.ii = [0,10];
% Linked Files and Directories
runopts.lflist = {fnameData};

%% Constructing bgs

% Bond Connectivity Map
bgs.bonds = uint16([3,1;
                    6,3;
                    4,3;
                    7,4;
                    2,7;
                    7,5;
                    8,7;
                    10,8;
                    10,9;
                    11,10]);

% Element List
bgs.elements = preallocEL(11);

bgs.elements(1).etype = 'MSF';
bgs.elements(1).exprtype = 'NMC';
bgs.elements(1).mvars = sym('t');
bgs.elements(1).expr = Qf;
                        
bgs.elements(2).etype = 'C';
bgs.elements(2).exprtype = 'CC';
bgs.elements(2).expr = C;
bgs.elements(2).x0 = 0.5;

bgs.elements(3).etype = '0';

bgs.elements(4).etype = 'R2';
bgs.elements(4).exprtype = 'NE';
bgs.elements(4).caus = int8([1,1]);
bgs.elements(4).expr = {Q1,Q1};

bgs.elements(5) = bgs.elements(1);

bgs.elements(6) = bgs.elements(2);
bgs.elements(6).x0 = 0.3;

bgs.elements(7) = bgs.elements(3);

bgs.elements(8).etype = 'MRN';
bgs.elements(8).exprtype = 'NMC';
bgs.elements(8).epnum = uint16(2);
bgs.elements(8).ebcon = uint16([7,8]);
bgs.elements(8).caus = int8([1,1]);
bgs.elements(8).mvars = sym('e6');
bgs.elements(8).expr = {Q2cond,Q2condMin;
                        Q2cond,Q2condMin};

bgs.elements(9) = bgs.elements(1);

bgs.elements(10) = bgs.elements(3);

bgs.elements(11) = bgs.elements(6);
bgs.elements(11).x0 = 0.1;

%% Storing runopts and bgs

save(fnameBGSD,'runopts','bgs');

%% Clearing memory

clear('C','Qa','omega','R1','R2','R2gr','Qf','Q1','Q2cond','Q2condMin',...
   'fnameData','fnameBGSD','runopts','bgs');