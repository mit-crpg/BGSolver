%% createxdotsys
% State derivative vector system creation function.
%
% xdotsys = createxdotsys(srtsys,vars,runopts)
%
% This function creates the state derivative vector system data structure, which
% is used by the function evalxdot to evaluate the state derivative vector. The
% way the state derivative functions are implemented depends on the options
% supplied in runtime and recording options data structure.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Variable descriptions
%
% BVID/NID means Bond Variable ID or Numeric Variable ID, which refers to
% absolute IDs of corresponding variable types. BVIDRFX/NIDRFX means Bond
% Variable ID Required For X-dot / Numeric Variable ID Required For X-dot, which
% refers to IDs of corresponding variables within the variable vectors that
% contain only the BVs/NVs required for x-dot.
% All variables that refer to these Required for X-dot variables have 'rfx'
% suffixes in their names.
%
% Numeric variables are separated into 7 types, based on their input types:
%        If reg. input        If time mod.         If mod. input
% Type 1      Yes                  No                   No
% Type 2      No                   Yes                  No
% Type 3      No                   No                   Yes
% Type 4      No                   Yes                  Yes
% Type 5      Yes                  No                   Yes
% Type 6      Yes                  Yes                  No
% Type 7      Yes                  Yes                  Yes
%
% During evaluation, the evaluation function loops over each NV type
% individually, because speed tests have shown this approach to be slightly
% faster than looping through the full NVRFX vector and checking NV type for
% every NV.

%% State derivative vector system creation function

function xdotsys = createxdotsys(srtsys,vars,runopts,ti_info)
% The purpose of this function is to create the state derivative vector system
% data structure, based on the sorting system and the runtime and recording
% options.
%
% INPUTS/OUTPUTS
% srtsys - Sorting system data structure.
% vars - Variables data structure.
% runopts - Runtime and Recording Options data structure.
% xdotsys - State derivative vector system data structure.
% ti_info - Time Integrator database data structure. Loaded from ti_info.mat.
%           OPTIONAL.

% Identifying if evaluator vectorization is required
if nargin == 4 && isfield(ti_info,runopts.tit)
   ti = ti_info.(runopts.tit);
   xdotsys.ifvectreq = ti.ifvectreq;
else
   xdotsys.ifvectreq = true;
end;

% Sizing the problem
bvnum = uint16(length(srtsys.b));
nvnum = uint16(length(srtsys.nhs));
xdotsys.xnum = srtsys.xnum;
xdotsys.bvnumrfx = uint16(length(srtsys.reqbvsxdot));
xdotsys.nvnumrfx = uint16(length(srtsys.reqnvsxdot));
xdotsys.lnumrfx = srtsys.lnumxdot;

% Identifying numeric variables and numeric function handles vectors
nvs = vars.n(srtsys.reqnvsxdot);
xdotsys.nhs = srtsys.nhs(srtsys.reqnvsxdot);

% Retrieving x-dot rfx BVIDs
xdotsys.xdotIDsrfx = getrfvids(srtsys.xdotIDs,srtsys.reqbvsxdot,bvnum);

% Preallocating layer-distributed rfx BVIDs and NIDs
xdotsys.lbvidsrfx = cell(1,srtsys.lnumxdot);
% Preallocating bv input variable ID arrays
xdotsys.lbvinpxids = cell(1,srtsys.lnumxdot); % XIDs input to BV layers
xdotsys.lbvinpnidsrfx = cell(1,srtsys.lnumxdot); % NIDs (RFX) input to BV layers
% Preallocating bv input statuses
xdotsys.bvifit = false(1,srtsys.lnumxdot); % If time is input to bv levels
xdotsys.bvifix = false(1,srtsys.lnumxdot); % If x is input to bv levels
xdotsys.bvifin = false(1,srtsys.lnumxdot); % If n is input to bv levels
% Preallocating bv in-input list positions
% Time position in bv input
xdotsys.bvtPosInInp = zeros(1,srtsys.lnumxdot,'uint8');
% State position in bv input
xdotsys.bvxPosInInp = zeros(1,srtsys.lnumxdot,'uint8');
% NV position in bv input
xdotsys.bvnPosInInp = zeros(1,srtsys.lnumxdot,'uint8');
% Preallocating bv input list size vector
xdotsys.bvInpSiz = zeros(1,srtsys.lnumxdot,'uint8');
% Preallocating nvidrfx array per layer and per nv type
xdotsys.lnvidsrfxByType = cell(7,srtsys.lnumxdot-1);
% Preallocating nv input tracking array per layer and per nv type
xdotsys.lnvinpsByType = cell(7,srtsys.lnumxdot-1);

% Combining regular and modulating input tracking arrays
nvinpids = [srtsys.nvinpxids(srtsys.reqnvsxdot,:) ...
   srtsys.nvinpbvids(srtsys.reqnvsxdot,srtsys.reqbvsxdot)];
nvminpids = [srtsys.nvminpxids(srtsys.reqnvsxdot,:) ...
	srtsys.nvminpbvids(srtsys.reqnvsxdot,srtsys.reqbvsxdot)];
% Filling nvrfx input statuses
nvifinp = full(any(nvinpids,2));
nvifmit = srtsys.nvmitids(srtsys.reqnvsxdot);
nvifminp = full(any(nvminpids,2));
% Filling nvrfx nv type statuses
nviftyp = cell(7,1);
nviftyp{1} = nvifinp & ~nvifmit & ~nvifminp;
nviftyp{2} = ~nvifinp & nvifmit & ~nvifminp;
nviftyp{3} = ~nvifinp & ~nvifmit & nvifminp;
nviftyp{4} = ~nvifinp & nvifmit & nvifminp;
nviftyp{5} = nvifinp & ~nvifmit & nvifminp;
nviftyp{6} = nvifinp & nvifmit & ~nvifminp;
nviftyp{7} = nvifinp & nvifmit & nvifminp;

% Converting the bv expressions to layer-distributed vector functions
xdotsys.b = cell(1,srtsys.lnumxdot); % BV (RFX) expression vector functions
for lid = 1:srtsys.lnumxdot
   % Retrieving layer's BVs and BVIDs vectors
   lbvids = srtsys.lbvsrfxdot{lid}; % Absolute BVIDs
   xdotsys.lbvidsrfx{lid} = getrfvids(lbvids,srtsys.reqbvsxdot,bvnum);
   lbvs = srtsys.b(lbvids);
   % Retrieving layer's bv variable dependencies tracking arrays
   lbvitids = find(srtsys.bvitids(lbvids));
   xdotsys.lbvinpxids{lid} = uint16(find(any(srtsys.bvinpxids(lbvids,:),1)));
   xdotsys.lbvinpnidsrfx{lid} = ...
      uint16(find(any(srtsys.bvinpnids(lbvids,srtsys.reqnvsxdot),1)));
   % Filling layer's bv input statuses
   xdotsys.bvifit(lid) = any(lbvitids);
   xdotsys.bvifix(lid) = any(xdotsys.lbvinpxids{lid});
   xdotsys.bvifin(lid) = any(xdotsys.lbvinpnidsrfx{lid});
   % Filling layer's bv in-input array positions
   xdotsys.bvtPosInInp(lid) = xdotsys.bvifit(lid);
   xdotsys.bvxPosInInp(lid) = xdotsys.bvifix(lid) * (xdotsys.bvifit(lid) + ...
      xdotsys.bvifix(lid));
   xdotsys.bvnPosInInp(lid) = xdotsys.bvifin(lid) * (xdotsys.bvifit(lid) + ...
      xdotsys.bvifix(lid) + xdotsys.bvifin(lid));
   % Computing number of bv in-input array positions
   xdotsys.bvInpSiz(lid) = ...
      xdotsys.bvifit(lid) + xdotsys.bvifix(lid) + xdotsys.bvifin(lid);
   
   % Checking for potential scalar BVs
   if lid == 1
      scalarLBVIDs = uint16(find(~(full(any(srtsys.bvinpxids(lbvids,:),2)) | ...
         full(any(srtsys.bvinpnids(lbvids,srtsys.reqnvsxdot),2)))));
      if isempty(scalarLBVIDs)
         xdotsys.bv1vectstat = int8(0); % No vectorization issues
      elseif length(scalarLBVIDs) < length(xdotsys.lbvidsrfx{1})
         xdotsys.bv1vectstat = int8(1); % .m function vectorization required
      elseif length(scalarLBVIDs) == length(xdotsys.lbvidsrfx{1})
         xdotsys.bv1vectstat = int8(2); % BV layer 1 fully scalar
      end;
   end;
   
   % Checking if symbolic expressions must be converted to numeric
   if runopts.cs2n
      % Checking numeric function type
      if runopts.nft == 1 % Temporary MATLAB .m function files
         % Setting temporary file path
         fpath = sprintf('TmpFunctions/bvrfx%u.m',lid);
         % Constructing input variables array
         invars = cell(1,xdotsys.bvInpSiz(lid));
         if xdotsys.bvifit(lid)
            invars{xdotsys.bvtPosInInp(lid)} = vars.ui.t;
         end;
         if xdotsys.bvifix(lid)
            invars{xdotsys.bvxPosInInp(lid)} = vars.x(xdotsys.lbvinpxids{lid});
         end;
         if xdotsys.bvifin(lid)
            invars{xdotsys.bvnPosInInp(lid)} = nvs(xdotsys.lbvinpnidsrfx{lid});
         end;
         % Constructing temporary MATLAB .m function
         xdotsys.b{lid} = matlabFunction(lbvs,'file',fpath,'vars',invars);
         % Vectorizing .m function, if necessary
         if lid == 1 && xdotsys.bv1vectstat == 1
            vectorizeBVfunction(fpath,scalarLBVIDs);
         end;
      else
         error('Only runopts.nft = 1 is presently supported!');
      end;
   else
      error('In this version, runopts.cs2n must be on!');
   end;
   
   % Checking for top layer (no numeric variables possible)
   if lid < srtsys.lnumxdot
      % Retrieving layer's NVID RFX vector
      lnids = srtsys.lnvsrfxdot{lid};
      [~,lnidsrfxBool] = getrfvids(lnids,srtsys.reqnvsxdot,nvnum,'vert');
      
      % Identifying layer's NVIDS RFX of all types
      for nvtyp = 1:7
         xdotsys.lnvidsrfxByType{nvtyp,lid} = ...
            uint16(find(lnidsrfxBool & nviftyp{nvtyp})).';
      end;
      
      % Filling out NV RFX input tracking arrays
      % NV type 1
      xdotsys.lnvinpsByType{1,lid}.nvinpids = ...
         getrightindvects(nvinpids,xdotsys.lnvidsrfxByType{1,lid});
      % NV type 2
      % Type 2 NVs are only time-modulated
      % NV type 3
      xdotsys.lnvinpsByType{3,lid}.nvminpids = ...
         getrightindvects(nvminpids,xdotsys.lnvidsrfxByType{3,lid});
      % NV type 4
      xdotsys.lnvinpsByType{4,lid}.nvminpids = ...
         getrightindvects(nvminpids,xdotsys.lnvidsrfxByType{4,lid});
      % NV type 5
      xdotsys.lnvinpsByType{5,lid}.nvinpids = ...
         getrightindvects(nvinpids,xdotsys.lnvidsrfxByType{5,lid});
      xdotsys.lnvinpsByType{5,lid}.nvminpids = ...
         getrightindvects(nvminpids,xdotsys.lnvidsrfxByType{5,lid});
      % NV type 6
      xdotsys.lnvinpsByType{6,lid}.nvinpids = ...
         getrightindvects(nvinpids,xdotsys.lnvidsrfxByType{6,lid});
      % NV type 7
      xdotsys.lnvinpsByType{7,lid}.nvinpids = ...
         getrightindvects(nvinpids,xdotsys.lnvidsrfxByType{7,lid});
      xdotsys.lnvinpsByType{7,lid}.nvminpids = ...
         getrightindvects(nvminpids,xdotsys.lnvidsrfxByType{7,lid});
   end;
end;

% Rehashing MATLAB path
rehash();
end