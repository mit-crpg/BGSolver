%% BGSD_defaultsBuild
% .mat file building script for the BGSD_Generator subpackage.
% 
% This script builds the defaults .mat file that is used by the
% BGSD_Generator subpackage, and stores it in the mat files folder. This
% .mat file is called 'BGSD_defaults.mat'.
%
% Package:    BGSolver v1.03
% Subpackage: mat_Files
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Filling runopts data structure

% SECTION 1 - File Header.
defaults.runopts.svc = 1.03;
defaults.runopts.svf = 1.03;
defaults.runopts.ifnotes = false;
defaults.runopts.notes = 'Default notes text.';

% SECTION 2 - Runtime and Recording Options.
defaults.runopts.verbose = int8(1);
defaults.runopts.rti = true;
defaults.runopts.log = int8(0);
defaults.runopts.logfnspec = int8(0);
defaults.runopts.smr = true;
defaults.runopts.smfnspec = int8(0);
defaults.runopts.rtrr = false;
defaults.runopts.rtrrfqspec = int8(2);
defaults.runopts.rtrrfnspec = int8(0);
defaults.runopts.rrm = true;
defaults.runopts.rrmfnspec = int8(0);

% SECTION 3 - Time Integrator Options.
defaults.runopts.titspec = int8(1);
defaults.runopts.tit = 'ode15s';
defaults.runopts.iispec = int8(1);
defaults.runopts.tst = int8(0);
defaults.runopts.ftsspect = int8(0);
defaults.runopts.ftsspec = int8(1);
defaults.runopts.phxsplit = int8(0);
defaults.runopts.tidisp = int8(0);

% SECTION 4 - Solver and Jumpstart Options.
defaults.runopts.stolx = 0;
defaults.runopts.stolres = 0;
defaults.runopts.sdisp = int8(-1);
defaults.runopts.srscal = int8(-1);
defaults.runopts.sxscal = int8(-1);
defaults.runopts.sxcent = int8(-1);
defaults.runopts.spc = int8(-1);
defaults.runopts.sigpdct = int8(-1);
defaults.runopts.sjacob = int8(-1);
defaults.runopts.sjsp = int8(-1);
defaults.runopts.jstype = int8(-1);
defaults.runopts.jsmsnum = uint32(0);
defaults.runopts.jsphxsplit = int8(-1);
defaults.runopts.jstolx = 0;
defaults.runopts.jstolres = 0;
defaults.runopts.jstit = 'default';
defaults.runopts.jsdisp = int8(-1);
defaults.runopts.jsrscal = int8(-1);
defaults.runopts.jsxscal = int8(-1);
defaults.runopts.jsxcent = int8(-1);
defaults.runopts.jspc = int8(-1);
defaults.runopts.jsjacob = int8(-1);
defaults.runopts.jsjsp = int8(-1);

% SECTION 5 - Evaluator Options.
defaults.runopts.cs2n = true;
defaults.runopts.nft = int8(1);

% SECTION 6 - Post-Processor Options.
defaults.runopts.epp = int32(-1);
defaults.runopts.fpp = int32(-1);
defaults.runopts.ppdisp = int8(0);

% SECTION 7 - Linked Files and Directories.
defaults.runopts.lfnum = uint8(0);
defaults.runopts.ldnum = uint8(0);

%% Filling bgs data structure

% SECTION 8 - Additional System Information.
defaults.bgs.asi.sci = 'FC';
defaults.bgs.asi.smi = 'UMP';
defaults.bgs.asi.sli = 'UL';
defaults.bgs.asi.sri = 'UR';
defaults.bgs.asi.sei = 'DNCE';

% SECTION 9 - Bond Connectivity Map.
% ...

% SECTION 10 - Element List.
defaults.bgs.element.phx = uint8(0);
defaults.bgs.element.nn = false;
defaults.bgs.element.at = 1.0e-6;

%% Saving the data structures in a .mat file

load('BGSolver_paths.mat'); % Provides paths to subpackages

fname = 'BGSD_defaults.mat'; 
dirname = BGSolverPaths.mat_Files;
filename = fullfile(dirname,fname);
save(filename,'defaults');

%% Clearing memory

clear('defaults','BGSolverPaths','fname','dirname','filename');