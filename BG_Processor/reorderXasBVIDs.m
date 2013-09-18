%% reorderXasBVIDs
% State derivative vector reordering in order of increasing BVIDs function.
%
% [srtsys,vars,refs] = reorderXasBVIDs(srtsys,vars,refs)
%
% This function reorders the state vector by reordering the state derivative
% vector in order of increasing BVIDs, instead of in order of increasing EIDs of
% the storage elements associated with the state variables. It then propagates
% the reorder to the appropriate XID tracking and referencing arrays. The
% function assumes the XID tracking arrays in the sorting system to already be
% constructed and be in sparse matrix format.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Main state derivative vector reordering function code

function [srtsys,vars,refs] = reorderXasBVIDs(srtsys,vars,refs)
% The purpose of this function is to reorder the state vector to get the state
% derivative vector in order of increasing corresponding BVIDs, instead of in
% order of increasing EIDs of the storage elements associated with the state
% variables. The reordering is then propagated to appropriate XID tracking
% arrays.
%
% INPUTS/OUTPUTS
% srtsys - Sorting System data structure.
% vars - Variables data structure.
% refs - Referencing arrays data structure.
% See the variable description cell in fBGSolve for more information.

% Reordering the state derivative vector
[srtsys.xdotIDs,srtsys.newxinds] = sort(srtsys.xdotIDs);

% Propagating the reordering to vars
vars.x = vars.x(srtsys.newxinds);

% Propagating the reordering to refs
[~,oldrefxinds,oldrefxids] = find(refs.xids);
refs.xids(oldrefxinds) = oldrefxids(srtsys.newxinds);

% Propagating the reordering to srtsys
srtsys.x0s = srtsys.x0s(srtsys.newxinds);
srtsys.xphxs = srtsys.xphxs(srtsys.newxinds);
srtsys.xnns = srtsys.xnns(srtsys.newxinds);
srtsys.xats = srtsys.xats(srtsys.newxinds);
srtsys.nvinpxids = srtsys.nvinpxids(:,srtsys.newxinds);
srtsys.nvminpxids = srtsys.nvminpxids(:,srtsys.newxinds);
srtsys.nvdepxids = srtsys.nvdepxids(:,srtsys.newxinds);

% Computing reverse ordering vectors
[~,srtsys.revxinds] = sort(srtsys.newxinds);
srtsys.revxinds = uint16(srtsys.revxinds);
end