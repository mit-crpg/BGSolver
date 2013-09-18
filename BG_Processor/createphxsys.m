%% createphxsys
% Physics splitting system creation function.
%
% phxsys = createphxsys(srtsys,runopts)
%
% This function creates the physics splitting system data structure, which is
% used by the time integrator for operator splitting.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Physics splitting system creation function

function phxsys = createphxsys(srtsys,runopts)
% The purpose of this function is to create the physics splitting data
% structure, based on the sorting system and the runtime and recording options.
%
% INPUTS/OUTPUTS
% srtsys - Sorting system data structure.
% runopts - Runtime and Recording Options data structure.
% phxsys - Physics splitting system data structure.

% Sizing the problem
phxsys.phxnum = max(srtsys.xphxs);
phxsys.xnum = srtsys.xnum;
phxsys.xphxs = cell(1,phxsys.phxnum);

% Checking if physics are split
if phxsys.phxnum == 0
   phxsys.ifphxsplit = false;
else
   phxsys.ifphxsplit = true;
   xphxs = srtsys.xphxs;
   xphxs(xphxs == 0) = phxsys.phxnum;
   % Looping through physics levels
   for phxid = 1:phxsys.phxnum
      phxsys.xphxs{phxid} = uint16(find(xphxs == phxid));
   end;
end;
phxsys.phxsplit = runopts.phxsplit;
end