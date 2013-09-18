%% installBGSolver
% Installation script for BGSolver v1.03.
%
% This script assumes that the current directory is the BGSolver v1.03 Installer
% directory.
%
% This script does NOT request any locations from the user. It adds the
% locations of all subpackages of BGSolver v1.03 to the path. The locations are
% also saved in a .mat file which can be used by the uninstaller. The .mat file
% is placed in the current directory. If previous versions of the same
% subpackages are present on the path, this script does not remove them, and it
% is up to the user to manually clean the path of duplicates.
% This script does NOT compile any of the .c files in BGSolver v1.03; this is
% done by compileBGSolver script. This script DOES run the mat_Files
% construction scripts.
% 
% The following subpackages are installed:
% - Installer
% - BGSD_Generator
% - mat_Files
% - BG_Processor
% - Time_Integrators
% - Utilities
%
% Package:    BGSolver v1.03
% Subpackage: Installer
% Date:       April 18, 2013
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Inducing BGSolver's location

BGSolverPaths.root = fileparts(cd);

%% Constructing subpackages' paths

BGSolverPaths.Installer = [BGSolverPaths.root filesep 'Installer'];
BGSolverPaths.BGSD_Generator = [BGSolverPaths.root filesep ...
   'BGSD_Generator'];
BGSolverPaths.mat_Files = [BGSolverPaths.root filesep 'mat_Files'];
BGSolverPaths.mat_scripts = [BGSolverPaths.mat_Files filesep 'scripts'];
BGSolverPaths.BG_Processor = [BGSolverPaths.root filesep 'BG_Processor'];
BGSolverPaths.Time_Integrators = [BGSolverPaths.root filesep ...
   'Time_Integrators'];
BGSolverPaths.Utilities = [BGSolverPaths.root filesep 'Utilities'];
BGSolverPaths.GetFullPath = [BGSolverPaths.Utilities filesep 'GetFullPath'];
BGSolverPaths.MEXUtilities = [BGSolverPaths.Utilities filesep 'MEX Utilities'];

utilNames = {'Array Manipulation','Data Manipulation','File Operations',...
   'GetFullPath','MEX Utilities','String Utilities'};
utilPathFun = @(utilName) [BGSolverPaths.Utilities filesep utilName];
BGSolverPaths.SubUtilities = cellfun(utilPathFun,utilNames,...
   'UniformOutput',false);

%% Recording the locations in a .mat file

fname = 'BGSolver_paths.mat';
save(fname,'BGSolverPaths');

%% Adding the subpackages' paths to the search path

addpath(BGSolverPaths.Installer,BGSolverPaths.BGSD_Generator,...
   BGSolverPaths.mat_Files,BGSolverPaths.mat_scripts,...
   BGSolverPaths.BG_Processor,BGSolverPaths.Time_Integrators,...
   BGSolverPaths.Utilities,BGSolverPaths.SubUtilities{:},'-begin');
savepath;

%% Building mat_Files

runAllmatFileScripts;

%% Reporting installed subpackages

fprintf('BGSolver v1.03 package has been installed.\n');
fprintf('The following subpackages have been installed:\n');
fprintf('- Installer\n');
fprintf('- BGSD_Generator\n');
fprintf('- mat_Files\n');
fprintf('- BG_Processor\n');
fprintf('- Time_Integrators\n');
fprintf('- Utilities\n');

%% Clearing memory

clear('BGSolverPaths','fname','utilNames','utilPathFun');