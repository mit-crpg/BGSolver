%% getvars2evalpp
% Identification function for bond and numeric variables to evaluate for
% post-processing.
%
% srtsys = getvars2evalpp(srtsys,runopts)
%
% This function identifies which bond and numeric variables have to be evaluated
% for bond variable evaluations during post-processing.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master variable identification for evaluation function

function srtsys = getvars2evalpp(srtsys,runopts,refs)
% The purpose of this function is to identify which bond and numeric variables
% have to be evaluated for bond variable evaluations during post-processing.
%
% INPUTS/OUTPUTS
% srtsys - Sorting system data structure.
% runopts - Runtime and Recording Options data structure, which contains
%           information about post-processing.
% refs - Referencing arrays data structure.

% Constructing BV numeric layer tracking arrays
if runopts.epp == 0 && runopts.fpp == 0
   srtsys.ifpp = false;
else
   srtsys.ifpp = true;
   if runopts.epp == 0
      eppIDs = false(1,srtsys.bvnum);
      eppIDsOrdered = zeros(1,0,'uint16');
   elseif runopts.epp == -1
      eppIDs = false(1,srtsys.bvnum);
      eppIDs(refs.ebvid(1:srtsys.bnum)) = true;
      eppIDsOrdered = refs.ebvid(1:srtsys.bnum);
   elseif runopts.epp > 0
      eppIDs = false(1,srtsys.bvnum);
      eppIDs(refs.ebvid(runopts.eppIDs)) = true;
      eppIDsOrdered = refs.ebvid(runopts.eppIDs);
   end;
   if runopts.fpp == 0
      fppIDs = false(1,srtsys.bvnum);
      fppIDsOrdered = zeros(1,0,'uint16');
   elseif runopts.fpp == -1
      fppIDs = false(1,srtsys.bvnum);
      fppIDs(refs.fbvid(1:srtsys.bnum)) = true;
      fppIDsOrdered = refs.fbvid(1:srtsys.bnum);
   elseif runopts.fpp > 0
      fppIDs = false(1,srtsys.bvnum);
      fppIDs(refs.fbvid(runopts.fppIDs)) = true;
      fppIDsOrdered = refs.fbvid(runopts.fppIDs);
   end;
   ppIDs = eppIDs | fppIDs;
   srtsys.ppIDs = uint16(find(ppIDs));
   srtsys.lnumpp = max(srtsys.bvlayers(ppIDs));
   srtsys.layerbvspp = cell(srtsys.lnumpp,1);
end;
if srtsys.ifpp
   for lid = 1:srtsys.lnumpp
      bvsl = (srtsys.bvlayers == lid & ppIDs);
      srtsys.layerbvspp{lid} = uint16(find(bvsl));
   end;
end;

% Checking if any post-processing occurs
if srtsys.ifpp
   reqbvspp = false(1,srtsys.bvnum);
   reqbvspp(srtsys.layerbvspp{end}) = true;
   srtsys.lbvsrfpp = cell(srtsys.lnumpp,1);
   srtsys.lbvsrfpp{end} = srtsys.layerbvspp{end};
   reqnvspp = false(1,srtsys.nvnum);
   srtsys.lnvsrfpp = cell(srtsys.lnumpp-1,1);
   for lid = (srtsys.lnumpp-1):(-1):1
      % Identifying all NVs on the current layer
      nvsl = srtsys.layernvs{lid}; % NVs on the current layer
      % Identifying NVs required by the upper layers' BVs
      reqnvs = full(any(srtsys.bvinpnids(reqbvspp,:),1));
      % Selectring current layer's NVs required by the upper layers' BVs
      reqnvsl = false(1,srtsys.nvnum);
      reqnvsl(nvsl) = reqnvs(nvsl);
      % Converting logical indexing to integer indexing
      srtsys.lnvsrfpp{lid} = uint16(find(reqnvsl));
      % Combining required NVs with this layer's NVs
      reqnvspp = reqnvspp | reqnvsl;
      
      % Identifying all BVs on the current layer
      bvsl = srtsys.layerbvs{lid}; % BVs on the current layer
      % Identifying BVs required by the upper layers' NVs
      reqbvs = full(any(srtsys.nvdepbvids(reqnvspp,:),1));
      % Selecting current layer's BVs required by the upper layers' NVs
      reqbvsl = false(1,srtsys.bvnum);
      reqbvsl(bvsl) = reqbvs(bvsl);
      % Appending pp BVs on the current layer to the required BVs
      reqbvsl(srtsys.layerbvspp{lid}) = true;
      % Converting logical indexing to integer indexing
      srtsys.lbvsrfpp{lid} = uint16(find(reqbvsl));
      % Combining required BVs with this layer's BVs
      reqbvspp = reqbvspp | reqbvsl;
   end;
   srtsys.reqbvspp = uint16(find(reqbvspp));
   srtsys.reqnvspp = uint16(find(reqnvspp));
   srtsys.eppIDsRFPP = getrfvids(eppIDsOrdered,srtsys.reqbvspp,srtsys.bvnum);
   srtsys.fppIDsRFPP = getrfvids(fppIDsOrdered,srtsys.reqbvspp,srtsys.bvnum);
end;
end