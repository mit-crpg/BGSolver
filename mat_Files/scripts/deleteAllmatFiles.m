%% deleteAllmatFiles
% Script for deleting all .mat files in the mat_Files subpackage.
%
% This script deletes all .mat files in the mat_Files subpackage. It then
% provides output for the user about which files it could find and delete.
%
% Package:    BGSolver v1.03
% Subpackage: mat_Files
% Date:       April 18, 2013
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Loading file paths

load('BGSolver_paths.mat');
mfdirname = BGSolverPaths.mat_Files;

%% Setting file names

matfnames = {'BGSD_defaults.mat','ti_info.mat'};

%% Deleting .mat files individually

mdelstats = false(size(matfnames));

% Looping through the .mat files
for mfid = 1:length(matfnames)
   mfname = matfnames{mfid};
   mfilename = fullfile(mfdirname,mfname);
   if exist(mfilename,'file') == 2
      delete(mfilename);
      mdelstats(mfid) = true;
   end;
end;

%% Reporting deleted .mat files

if any(mdelstats)
   fprintf('Deleted the following .mat file(s):\n');
   for mfid = 1:length(mdelstats)
      if mdelstats(mfid)
         fprintf('- %s\n',matfnames{mfid});
      end;
   end;
end;

if any(~mdelstats)
   fprintf('The following .mat file(s) have not been found and deleted:\n');
   for mfid = 1:length(mdelstats)
      if ~mdelstats(mfid)
         fprintf('- %s\n',matfnames{mfid});
      end;
   end;
end;

%% Clearing memory

clear('BGSolverPaths','mfdirname','matfnames','mdelstats','mfid','mfname',...
   'mfilename');