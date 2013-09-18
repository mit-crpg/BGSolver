%% ti_infoBuild
% .mat file building script for the BG_Processor subpackage.
%
% This script builds the time integrator info .mat file that is used by the
% BG_Processor subpackage, and stores it in the mat files folder. This .mat file
% is called 'ti_info.mat'.
%
% Package:    BGSolver v1.03
% Subpackage: mat_Files
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% ode15s

ti_info.ode15s.src = 'Built-in MATLAB function';
ti_info.ode15s.stiff = 'Stiff';
ti_info.ode15s.order = 'Low to medium (varies)';
ti_info.ode15s.explicit = 'N/A';
ti_info.ode15s.method = 'Advanced';
ti_info.ode15s.scaling = 'Built-in MATLAB time integrator';
ti_info.ode15s.vector = 'Supported';
ti_info.ode15s.opsplit = 'Full Coupling only';
ti_info.ode15s.dt = 'Adaptive only';
ti_info.ode15s.rr = 'Not supported';
ti_info.ode15s.jstype = 'Jumpstart not needed';
ti_info.ode15s.rscal = 'Not supported';
ti_info.ode15s.xscal = 'Not supported';
ti_info.ode15s.xcent = 'Unknown if supported';
ti_info.ode15s.pc = 'Unknown if supported';
ti_info.ode15s.igpdct = 'Unknown if supported';
ti_info.ode15s.jacob = 'numjac-based, sparse and vectorized';
ti_info.ode15s.addinfo = 'Most versatile integrator';
ti_info.ode15s.ifopsplit = false;
ti_info.ode15s.ifvectreq = true;
ti_info.ode15s.Njs = 0;
ti_info.ode15s.pts = false;

%% ode45

ti_info.ode45.src = 'Built-in MATLAB function';
ti_info.ode45.stiff = 'Nonstiff';
ti_info.ode45.order = 'Medium (varies)';
ti_info.ode45.explicit = 'Explicit';
ti_info.ode45.method = 'Runge-Kutta';
ti_info.ode45.scaling = 'Built-in MATLAB time integrator';
ti_info.ode45.vector = 'Allowed, but not used';
ti_info.ode45.opsplit = 'Full Coupling only';
ti_info.ode45.dt = 'Adaptive only';
ti_info.ode45.rr = 'Not supported';
ti_info.ode45.jstype = 'Jumpstart not needed';
ti_info.ode45.rscal = 'No residual used (explicit time integrator)';
ti_info.ode45.xscal = 'Not applicable';
ti_info.ode45.xcent = 'Unknown if supported';
ti_info.ode45.pc = 'No nonlinear solves used (explicit time integrator)';
ti_info.ode45.igpdct = 'Unknown if supported';
ti_info.ode45.jacob = 'Not used (explicit time integrator)';
ti_info.ode45.addinfo = 'MATLAB''s fastest integrator for nonstiff problems';
ti_info.ode45.ifopsplit = false;
ti_info.ode45.ifvectreq = false;
ti_info.ode45.Njs = 0;
ti_info.ode45.pts = false;

% Insert more time integrators' information here

%% Saving the data structure in a .mat file

load('BGSolver_paths.mat'); % Provides paths to subpackages

fname = 'ti_info.mat';
dirname = BGSolverPaths.mat_Files;
filename = fullfile(dirname,fname);
save(filename,'ti_info');

%% Clearing memory

clear('ti_info','BGSolverPaths','fname','dirname','filename');