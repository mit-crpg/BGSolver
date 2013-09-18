%% createppsys
% Post-processing vector system creation function.
%
% ppsys = createppsys(srtsys,vars,runopts)
%
% This function creates the post-processing vector system data structure, which
% is used by the function evalpp to evaluate the post-processing bond variables.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Post-processing vector system creation function.

function ppsys = createppsys(srtsys,vars,runopts)
% The purpose of this function is to create the post-processing vector system
% data structure, based on the sorting system and the runtime and recording
% options.
%
% INPUTS/OUTPUTS
% srtsys - Sorting system data structure.
% vars - Variables data structure.
% runopts - Runtime and Recording Options data structure.
% ppsys - Post-processing vector system data structure.

ppsys.ifpp = srtsys.ifpp;

if ~ppsys.ifpp
   return;
end;

% Setting the post-processing display status
ppsys.ppdisp = runopts.ppdisp;

% Sizing the problem
bvnum = uint16(length(srtsys.b));
nvnum = uint16(length(srtsys.nhs));
ppsys.xnum = srtsys.xnum;
ppsys.bvnumrfpp = uint16(length(srtsys.reqbvspp));
ppsys.nvnumrfpp = uint16(length(srtsys.reqnvspp));
ppsys.lnumrfpp = srtsys.lnumpp;

% Identifying nv numeric function handles and variables vectors
nvs = vars.n(srtsys.reqnvspp);
ppsys.nhs = srtsys.nhs(srtsys.reqnvspp);

% Retrieving Effort and Flow BVIDs RFPP
ppsys.eppidsrfpp = srtsys.eppIDsRFPP;
ppsys.fppidsrfpp = srtsys.fppIDsRFPP;

% Preallocating layer-distributed rfpp BVIDs and NIDs
ppsys.lbvidsrfpp = cell(1,srtsys.lnumpp);
% Preallocating bv input variable ID arrays
ppsys.lbvinpxids = cell(1,srtsys.lnumpp); % XIDs input to BV layers
ppsys.lbvinpnidsrfpp = cell(1,srtsys.lnumpp); % NIDs (RFPP) input to BV layers
% Preallocating bv input statuses
ppsys.bvifit = false(1,srtsys.lnumpp); % If time is input to bv levels
ppsys.bvifix = false(1,srtsys.lnumpp); % If x is input to bv levels
ppsys.bvifin = false(1,srtsys.lnumpp); % If n is input to bv levels
% Preallocating bv in-input list positions
ppsys.bvtPosInInp = zeros(1,srtsys.lnumpp,'uint8'); % Time position in bv input
ppsys.bvxPosInInp = zeros(1,srtsys.lnumpp,'uint8'); % State position in bv input
ppsys.bvnPosInInp = zeros(1,srtsys.lnumpp,'uint8'); % NV position in bv input
% Preallocating bv input list size vector
ppsys.bvInpSiz = zeros(1,srtsys.lnumpp,'uint8');
% Preallocating nvidrfpp array per layer and per nv type
ppsys.lnvidsrfppByType = cell(7,srtsys.lnumpp-1);
% Preallocating nv input tracking array per layer and per nv type
ppsys.lnvinpsByType = cell(7,srtsys.lnumpp-1);

% Combining regular and modulating input tracking arrays
nvinpids = [srtsys.nvinpxids(srtsys.reqnvspp,:) ...
   srtsys.nvinpbvids(srtsys.reqnvspp,srtsys.reqbvspp)];
nvminpids = [srtsys.nvminpxids(srtsys.reqnvspp,:) ...
	srtsys.nvminpbvids(srtsys.reqnvspp,srtsys.reqbvspp)];
% Filling nvrfpp input statuses
nvifinp = full(any(nvinpids,2));
nvifmit = srtsys.nvmitids(srtsys.reqnvspp);
nvifminp = full(any(nvminpids,2));
% Filling nvrfpp nv type statuses
nviftyp = cell(7,1);
nviftyp{1} = nvifinp & ~nvifmit & ~nvifminp;
nviftyp{2} = ~nvifinp & nvifmit & ~nvifminp;
nviftyp{3} = ~nvifinp & ~nvifmit & nvifminp;
nviftyp{4} = ~nvifinp & nvifmit & nvifminp;
nviftyp{5} = nvifinp & ~nvifmit & nvifminp;
nviftyp{6} = nvifinp & nvifmit & ~nvifminp;
nviftyp{7} = nvifinp & nvifmit & nvifminp;

% Converting the bv expressions to layer-distributed vector functions
ppsys.b = cell(1,srtsys.lnumpp); % BV (RFPP) expression vector functions
for lid = 1:srtsys.lnumpp
   % Retrieving layer's BVs and BVIDs vectors
   lbvids = srtsys.lbvsrfpp{lid}; % Absolute BVIDs
   ppsys.lbvidsrfpp{lid} = getrfvids(lbvids,srtsys.reqbvspp,bvnum);
   lbvs = srtsys.b(lbvids);
   % Retrieving layer's bv variable dependencies tracking arrays
   lbvitids = find(srtsys.bvitids(lbvids));
   ppsys.lbvinpxids{lid} = uint16(find(any(srtsys.bvinpxids(lbvids,:),1)));
   ppsys.lbvinpnidsrfpp{lid} = ...
      uint16(find(any(srtsys.bvinpnids(lbvids,srtsys.reqnvspp),1)));
   % Filling layer's bv input statuses
   ppsys.bvifit(lid) = any(lbvitids);
   ppsys.bvifix(lid) = any(ppsys.lbvinpxids{lid});
   ppsys.bvifin(lid) = any(ppsys.lbvinpnidsrfpp{lid});
   % Filling layer's bv in-input array positions
   ppsys.bvtPosInInp(lid) = ppsys.bvifit(lid);
   ppsys.bvxPosInInp(lid) = ppsys.bvifix(lid) * (ppsys.bvifit(lid) + ...
      ppsys.bvifix(lid));
   ppsys.bvnPosInInp(lid) = ppsys.bvifin(lid) * (ppsys.bvifit(lid) + ...
      ppsys.bvifix(lid) + ppsys.bvifin(lid));
   % Computing number of bv in-input array positions
   ppsys.bvInpSiz(lid) = ...
      ppsys.bvifit(lid) + ppsys.bvifix(lid) + ppsys.bvifin(lid);
   
   % Checking if symbolic expressions must be converted to numeric
   if runopts.cs2n
      % Checking numeric function type
      if runopts.nft == 1 % Temporary MATLAB .m function files
         % Setting temporary file path
         fpath = sprintf('TmpFunctions/bvrfpp%u.m',lid);
         % Constructing input variables array
         invars = cell(1,ppsys.bvInpSiz(lid));
         if ppsys.bvifit(lid)
            invars{ppsys.bvtPosInInp(lid)} = vars.ui.t;
         end;
         if ppsys.bvifix(lid)
            invars{ppsys.bvxPosInInp(lid)} = vars.x(ppsys.lbvinpxids{lid});
         end;
         if ppsys.bvifin(lid)
            invars{ppsys.bvnPosInInp(lid)} = nvs(ppsys.lbvinpnidsrfpp{lid});
         end;
         % Constructing temporary MATLAB .m function
         ppsys.b{lid} = matlabFunction(lbvs,'file',fpath,'vars',invars);
      else
         error('Only runopts.nft = 1 is presently supported!');
      end;
   else
      error('In this version, runopts.cs2n must be on!');
   end;
   
   % Checking for top layer (no numeric variables possible)
   if lid < srtsys.lnumpp
      % Retrieving layer's NVID RFPP vector
      lnids = srtsys.lnvsrfpp{lid};
      [~,lnidsrfppBool] = getrfvids(lnids,srtsys.reqnvspp,nvnum,'vert');
      
      % Identifying layer's NVIDS RFPP of all types
      for nvtyp = 1:7
         ppsys.lnvidsrfppByType{nvtyp,lid} = ...
            uint16(find(lnidsrfppBool & nviftyp{nvtyp})).';
      end;
      
      % Filling out NV RFPP input tracking arrays
      % NV type 1
      ppsys.lnvinpsByType{1,lid}.nvinpids = ...
         getrightindvects(nvinpids,ppsys.lnvidsrfppByType{1,lid});
      % NV type 2
      % Type 2 NVs are only time-modulated
      % NV type 3
      ppsys.lnvinpsByType{3,lid}.nvminpids = ...
         getrightindvects(nvminpids,ppsys.lnvidsrfppByType{3,lid});
      % NV type 4
      ppsys.lnvinpsByType{4,lid}.nvminpids = ...
         getrightindvects(nvminpids,ppsys.lnvidsrfppByType{4,lid});
      % NV type 5
      ppsys.lnvinpsByType{5,lid}.nvinpids = ...
         getrightindvects(nvinpids,ppsys.lnvidsrfppByType{5,lid});
      ppsys.lnvinpsByType{5,lid}.nvminpids = ...
         getrightindvects(nvminpids,ppsys.lnvidsrfppByType{5,lid});
      % NV type 6
      ppsys.lnvinpsByType{6,lid}.nvinpids = ...
         getrightindvects(nvinpids,ppsys.lnvidsrfppByType{6,lid});
      % NV type 7
      ppsys.lnvinpsByType{7,lid}.nvinpids = ...
         getrightindvects(nvinpids,ppsys.lnvidsrfppByType{7,lid});
      ppsys.lnvinpsByType{7,lid}.nvminpids = ...
         getrightindvects(nvminpids,ppsys.lnvidsrfppByType{7,lid});
   end;
end;

% Rehashing MATLAB path
rehash();
end