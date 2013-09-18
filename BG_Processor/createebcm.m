%% createebcm
% Element-indexed bond connectivity map creation function.
%
% ebcm = createebcm(bgs)
%
% This function creates an element-indexed bond connectivity map data structure,
% for easily identifying bonds connected to a given element, or pointing to or
% from it.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master element-indexed bond connectivity map creation function.

function ebcm = createebcm(bgs)
% The purpose of this function is to create an element-indexed bond connectivity
% map data structure.
%
% INPUTS/OUTPUTS
% bgs - Bond Graph System data structure. Assumes that it has been completed.
% ebcm - Element-indexed Bond Connectivity Map data structure. Its fields are:
%  .to{EID} - Vector of BIDs of bonds pointing to element EID. uint16 vector.
%  .from{EID} - Vector of BIDs of bonds pointing to element EID. uint16 vector.
%  .conn{EID} - Vector of BIDs of bonds connected to element EID. uint16 vector.
% See the variable description cell in fBGSolve for more information.

% Preallocating arrays
ebcm.to = cell(bgs.enum,1);
ebcm.from = cell(bgs.enum,1);
ebcm.conn = cell(bgs.enum,1);

% Identifying connected bonds for each element
for eid = 1:bgs.enum
   ebcm.to{eid} = uint16(find(bgs.bonds(:,1) == eid)).';
   ebcm.from{eid} = uint16(find(bgs.bonds(:,2) == eid)).';
   ebcm.conn{eid} = sort([ebcm.to{eid},ebcm.from{eid}]);
end;
end