%% compileBGSolver
% Script for compiling the .c files in BGSolver v1.03 into MATLAB executables.
%
% All .c files in BGSolver v1.03 have associated .m function files; the package
% therefore can run without compiling any of the .c files. It is recommended to
% do so though, because .mex files generally perform significantly better.
% A configured MEX compiler is required prior to running this script - otherwise
% it is not possible to build any .mex files. See 'mex -setup' command for more
% information.
% All .mex files in BGSolver v1.03 were compiled for Windows x64 only (.mexw64
% extensions), and are therefore unfit for other platforms. The .c files can be
% compiled for those platforms though. Note, that OpenMP 2.0 support is
% necessary for compiling some of the .c files in BGSolver v1.03.
% This script assumes that BGSolver v1.03 has been installed, either manually,
% or through the installBGSolver script.
%
% Package:    BGSolver v1.03
% Subpackage: Installer
% Date:       April 15, 2013
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

%% .c file paths

cFPaths = {[BGSolverPaths.BG_Processor filesep 'evalxdot.c'],...
           [BGSolverPaths.GetFullPath filesep 'GetFullPath.c']};

% Relative .c file path construction function
chID = length(BGSolverPaths.root) + 1;
fhPath = @(cFPath) cFPath(chID:end);
relcFPaths = cellfun(fhPath,cFPaths,'UniformOutput',false);

%% Compiling .mex files

% Compilation status vector
mcompstats = true(1,length(cFPaths));
% MEX Utilities library header and source files
cfNameMexUtilities = [BGSolverPaths.MEXUtilities filesep 'mexUtilities.c'];

% Looping through .c files, compiling one at a time
for cfID = 1:length(cFPaths)
   mexCmdStr = 'mex -O -I"%s" "%s" "%s"';
   [cFPath,cFName,cFExt] = fileparts(cFPaths{cfID});
   
   fprintf('Compiling %s...\n',relcFPaths{cfID});
   oldpath = cd(cFPath);
   mexCmdStr = sprintf(mexCmdStr,BGSolverPaths.MEXUtilities,[cFName cFExt],...
      cfNameMexUtilities);
   try
      eval(mexCmdStr);
   catch err
      fprintf('Failed to compile %s.\n',relcFPaths{cfID});
      mcompstats(cfID) = false;
   end;
   cd(oldpath);
end;

%% Reporting compiled .mex files

if any(mcompstats)
   fprintf('The following files have been compiled into .%s:\n',mexext);
   fprintf('- %s\n',relcFPaths{mcompstats})
end;
if any(~mcompstats)
   fprintf('The following files have not been compiled into .%s:\n',mexext);
   fprintf('- %s\n',relcFPaths{~mcompstats});
end;

%% Clearing memory

clear('fname','BGSolverPaths','cFPaths','chID','fhPath','relcFPaths',...
   'mcompstats','cfNameMexUtilities','cfID','mexCmdStr','cFPath','cFName',...
   'cFExt','oldpath');