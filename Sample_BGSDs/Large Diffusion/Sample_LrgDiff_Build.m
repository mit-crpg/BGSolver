%% Sample_LrgDiff_Build 
% .mat file building script for the large diffusion BGSD in the Sample_BGSDs
% subpackage.
%
% This script builds the .mat file which contains runopts and bgs data
% structures for a linear primitive diffusion BGS with variable element number.
% The object of this test is to time the sorting step. All elements are CC.
%
% Package:    BGSolver v1.03
% Subpackage: Sample_BGSDs
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Setting parameters

% Material properties
rhoCp = 1;
k = 0.05;
% Geometry
L = 1; % Total domain length
% Initial temperature peak value (sinusoidal, with 1/2 period length of L,
% peaking in the middle of domain)
Tpk = 10;
% Boundary temperatures
Tbnd = [1,1];

% Discretization
N = 300; % Number of interior points
Xpts = linspace(0,L,N+2); % X-coordinates of all points

fnameData = 'LrgDiff_Data.mat';
fnameBGSD = 'LrgDiff_BGSD.mat';

%% Constructing runopts

% Runtime and Recording Options
runopts.verbose = 2;
runopts.log = 2;
% Time Integrator Options
runopts.ii = [0,10];
runopts.tidisp = 1;
% Linked Files and Directories
runopts.lflist = {fnameData};
% Post-processing (enabling only capacitors' efforts)
runopts.epp = N;
runopts.eppIDs = 5:4:((N-1)*4+5);
runopts.fpp = 0;

%% Constructing bgs

bnum = 4*(N+1)+1;
enum = (N+1)*4+2;

% Bond Connectivity Map
bgs.bonds = zeros(bnum,2,'uint16');
% Setting left BC
bgs.bonds(1,:) = [1,2]; % SE_0 -> 0_0
bgs.bonds(2,:) = [3,1]; % 0_0 -> 1_0,1
% Setting interior points
eid = 3;
for bid = 3:4:(bnum-3)
   bgs.bonds(bid,:) = [eid+1,eid]; % 1 -> R
   bgs.bonds(bid+1,:) = [eid+2,eid]; % 1 -> 0
   bgs.bonds(bid+2,:) = [eid+3,eid+2]; % 0 -> C
   bgs.bonds(bid+3,:) = [eid+4,eid+2]; % 0 -> 1
   eid = eid + 4;
end;
% Setting right BC
bgs.bonds(bnum-2,:) = [enum-2,enum-3]; % 1_N,N+1 -> R_N,N+1
bgs.bonds(bnum-1,:) = [enum-1,enum-3]; % 1_N,N+1 -> 0_N+1
bgs.bonds(bnum,:) = [enum,enum-1]; % 0_N+1 -> SE_N+1

% Element Expressions
h = L / (N+1);
R = h^2 / k;
C = 1 / rhoCp;

% Initial conditions
T0fun = @(x) (Tpk-Tbnd(1)) * sin(pi/L * x) + Tbnd(1);

% Element List
bgs.elements = preallocEL(enum);
% Setting left BC
bgs.elements(1).etype = '0';
bgs.elements(2).etype = 'SE';
bgs.elements(2).exprtype = 'CC';
bgs.elements(2).expr = Tbnd(1);
% Setting interior points
i = 2;
for eid = 3:4:(enum-4)
   bgs.elements(eid).etype = '1';
   
   bgs.elements(eid+1).etype = 'R';
   bgs.elements(eid+1).exprtype = 'CC';
   bgs.elements(eid+1).expr = R;
   
   bgs.elements(eid+2).etype = '0';
   
   bgs.elements(eid+3).etype = 'C';
   bgs.elements(eid+3).exprtype = 'CC';
   bgs.elements(eid+3).expr = C;
   bgs.elements(eid+3).x0 = T0fun(Xpts(i));
   i = i + 1;
end;
% Setting right BC
bgs.elements(enum-3).etype = '1';
bgs.elements(enum-2).etype = 'R';
bgs.elements(enum-2).exprtype = 'CC';
bgs.elements(enum-2).expr = R;
bgs.elements(enum-1).etype = '0';
bgs.elements(enum).etype = 'SE';
bgs.elements(enum).exprtype = 'CC';
bgs.elements(enum).expr = Tbnd(2);

%% Storing supporting quantities in a .mat file

save(fnameData,'rhoCp','k','L','Tpk','Tbnd','N','Xpts','h','R','C','T0fun');

%% Storing runopts and bgs

save(fnameBGSD,'runopts','bgs');

%% Clearing memory

clear('rhoCp','k','L','Tpk','Tbnd','N','Xpts','h','R','C','i','eid','bid',...
   'enum','bnum','T0fun','fnameData','fnameBGSD','runopts','bgs');