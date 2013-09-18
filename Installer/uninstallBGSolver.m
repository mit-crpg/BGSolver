%% uninstallBGSolver
% Uninstallation script for BGSolver v1.03.
%
% This script checks whether the BGSolver v1.03 package is installed. If it is,
% this script removes the locations of all of the BGSolver subpackages from the
% path, and deletes all the mat_Files that may be present.
%
% The following subpackages are uninstalled:
% - Installer
% - BGSD_Generator
% - mat_Files
% - BG_Processor
% - Time_Integrators
% - Utilities
%
% Package:    BGSolver v1.03
% Subpackage: Installer
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Retrieving the subpackages' locations

fname = 'BGSolver_paths.mat';
try
   load(fname);
catch err
   if strcmp(err.identifier,'MATLAB:load:couldNotReadFile')
      fprintf('An existing BGSolver v1.03 is not found.\n');
      clear('fname','err');
      return;
   else
      rethrow(err);
   end;
end;
dirname = BGSolverPaths.Installer;
filename = fullfile(dirname,fname);

%% Deleting individual .mat files from the mat_Files subpackage

deleteAllmatFiles;
load(fname); % Reloading BGSolverPaths

%% Removing the subpackages from the search path

delete(filename);
rmpath(BGSolverPaths.Installer,BGSolverPaths.BGSD_Generator,...
   BGSolverPaths.mat_Files,BGSolverPaths.mat_scripts,...
   BGSolverPaths.BG_Processor,BGSolverPaths.Time_Integrators,...
   BGSolverPaths.Utilities,BGSolverPaths.SubUtilities{:});
savepath;

%% Reporting uninstalled subpackages

fprintf('BGSolver v1.03 package has been uninstalled.\n');
fprintf('The following subpackages have been uninstalled:\n');
fprintf('- Installer\n');
fprintf('- BGSD_Generator\n');
fprintf('- mat_Files\n');
fprintf('- BG_Processor\n');
fprintf('- Time_Integrators\n');
fprintf('- Utilities\n');

%% Clearing memory

clear('BGSolverPaths','fname','dirname','filename');