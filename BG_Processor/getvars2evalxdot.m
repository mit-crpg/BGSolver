%% getvars2evalxdot
% Identification function for bond and numeric variables to evaluate for
% state derivative evaluation.
%
% srtsys = getvars2evalpp(srtsys)
%
% This function identifies which bond and numeric variables have to be evaluated
% for state derivative evaluation.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master variable identification for evaluation function

function srtsys = getvars2evalxdot(srtsys)
% The purpose of this function is to identify which bond and numeric variables
% have to be evaluated for state derivative evaluation.
%
% INPUTS/OUTPUTS
% srtsys - Sorting system data structure.

% Constructing BV numeric layer tracking arrays
srtsys.lnumxdot = max(srtsys.bvlayers(srtsys.xdotIDs));
xdotIDs = false(1,srtsys.bvnum);
xdotIDs(srtsys.xdotIDs) = true;
srtsys.layerbvsxdot = cell(srtsys.lnumxdot,1);
for lid = 1:srtsys.lnum
   bvsl = (srtsys.bvlayers == lid);
   if lid <= srtsys.lnumxdot
      srtsys.layerbvsxdot{lid} = uint16(find(bvsl & xdotIDs));
   end;
end;

% Identifying variables to evaluate for state derivative evaluation
reqbvsxdot = false(1,srtsys.bvnum);
reqbvsxdot(srtsys.layerbvsxdot{end}) = true;
srtsys.lbvsrfxdot = cell(srtsys.lnumxdot,1);
srtsys.lbvsrfxdot{end} = srtsys.layerbvsxdot{end};
reqnvsxdot = false(1,srtsys.nvnum);
srtsys.lnvsrfxdot = cell(srtsys.lnumxdot-1,1);
for lid = srtsys.lnumxdot-1:-1:1
   % Identifying all NVs on the current layer
   nvsl = srtsys.layernvs{lid}; % NVs on the current layer
   % Identifying NVs required by the upper layers' BVs
   reqnvs = full(any(srtsys.bvinpnids(reqbvsxdot,:),1));
   % Selecting current layer's NVs required by the upper layers' BVs
   reqnvsl = false(1,srtsys.nvnum);
   reqnvsl(nvsl) = reqnvs(nvsl);
   % Converting logical indexing to integer indexing
   srtsys.lnvsrfxdot{lid} = uint16(find(reqnvsl));
   % Combining required NVs with this layer's NVs
   reqnvsxdot = reqnvsxdot | reqnvsl;
   
   % Identifying all BVs on the current layer
   bvsl = srtsys.layerbvs{lid};  % BVs on the current layer
   % Identifying BVs required by the upper layers' NVs
   reqbvs = full(any(srtsys.nvdepbvids(reqnvsxdot,:),1));
   % Selecting current layer's BVs required by the upper layers' NVs
   reqbvsl = false(1,srtsys.bvnum);
   reqbvsl(bvsl) = reqbvs(bvsl);
   % Appending xdot BVs on the current layer to the required BVs
   reqbvsl(srtsys.layerbvsxdot{lid}) = true;
   % Converting logical indexing to integer indexing
   srtsys.lbvsrfxdot{lid} = uint16(find(reqbvsl));
   % Combining required BVs with this layer's BVs
   reqbvsxdot = reqbvsxdot | reqbvsl;
end;
srtsys.reqbvsxdot = uint16(find(reqbvsxdot));
srtsys.reqnvsxdot = uint16(find(reqnvsxdot));
end