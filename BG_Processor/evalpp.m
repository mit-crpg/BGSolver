%% evalpp
% Post-processing vectors evaluation function.
%
% [E,F] = evalpp(T,X,ppsys)
%
% This function evaluates the effort and flow vectors as instructed by
% post-processing inputs. The vectors are output as vertical (in time)
% horizontally concatenated vectors, ordered in order of increasing BID, left to
% right.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Post-processing vectors evaluation function

function [E,F] = evalpp(T,X,ppsys)
% The purpose of this function is to evaluate the effort and flow vectors as
% instructed by post-processing inputs. The vectors are output as vertical (in
% time) horizontally concatenated vectors, ordered in order of increasing BID,
% left to right.
%
% INPUTS/OUTPUTS
% T - Time vector. Vertical double vector.
% X - State matrix. Double matrix, vertical in time, horizontal in xid.
% ppsys - Post-processing vector system data structure.
% E - Effort matrix. Vertical in time, horizontal in BID. Only the efforts that
%     were requested for post-processing are output.
% F - Flow matrix. Vertical in time, horizontal in BID. Only the flows that were
%     requested for post-processing are output.

% Sizing the problem
tnum = length(T);

% Transposing the state vector
X = X.';

% Preallocating arrays
brfpp = zeros(ppsys.bvnumrfpp,tnum);
nrfpp = zeros(ppsys.nvnumrfpp,tnum);

% Looping through time steps
for tid = 1:tnum
   % Retrieving time and state
   t = T(tid);
   x = X(:,tid);
   
   % Going through every layer and evaluating bond and numeric variables
   for lid = 1:ppsys.lnumrfpp
      % Evaluating bond variables on the layer
      if ppsys.bvifit(lid) || ppsys.bvifix(lid) || ppsys.bvifin(lid)
         bvinvars = cell(1,ppsys.bvInpSiz(lid));
         if ppsys.bvifit(lid)
            bvinvars{ppsys.bvtPosInInp(lid)} = t;
         end;
         if ppsys.bvifix(lid)
            bvinvars{ppsys.bvxPosInInp(lid)} = x(ppsys.lbvinpxids{lid});
         end;
         if ppsys.bvifin(lid)
            bvinvars{ppsys.bvnPosInInp(lid)} = ...
               nrfpp(ppsys.lbvinpnidsrfpp{lid},tid);
         end;
         brfpp(ppsys.lbvidsrfpp{lid},tid) = ppsys.b{lid}(bvinvars{:});
      end;
      
      % Evaluating numeric variables on the layer
      if lid < ppsys.lnumrfpp
         % Constructing concatenated input vector
         c = [x;brfpp(:,tid)];
         % Looping through NV types
         for nvtyp = 1:7
            tnidsrfpp = ppsys.lnvidsrfppByType{nvtyp,lid};
            tnhs = ppsys.nhs(tnidsrfpp);
            tnvinps = ppsys.lnvinpsByType{nvtyp,lid};
            nrfpp(tnidsrfpp,tid) = evalnv(t,c,nvtyp,tnhs,tnvinps);
         end;
      end;
   end;
   
   % Outputting time step information
   if ppsys.ppdisp == 1
      ppstepdisp(t);
   end;
end;

% Transposing the BV RFPP matrix
brfpp = brfpp.';
% Constructing the effort and flow matrices
E = brfpp(:,ppsys.eppidsrfpp);
F = brfpp(:,ppsys.fppidsrfpp);
end

%% Numeric variable evaluation function

function n = evalnv(t,c,nvtyp,nhs,nvinps)
% The purpose of this function is to evaluate a vector of numeric variables of a
% given type.
%
% INPUTS/OUTPUTS
% t - Time to evaluate the numeric variables at.
% c - Concatenated input vector.
% nidsrfpp - NIDs RFPP of the numeric variables to evaluate.
% nvtype - Numeric variable type of the numeric variables to evaluate. See
%          createxdotsys variable documentation cell for details of NV types.
% nhs - Vector of numeric function handles to evaluate.
% nvinps - NV input tracking data structure.
% n - Evaluated numeric variable vector.

% Sizing the problem
nvnum = uint16(length(nhs));

% Prealloacting n
n = zeros(nvnum,1);

switch nvtyp
   % Type 1
   % - Regular input
   case 1
      % Looping through NVs to evaluate
      for nid = 1:nvnum
         n(nid) = nhs{nid}(c(nvinps.nvinpids{nid}));
      end;
   
   % Type 2
   % - Time modulation
   case 2
      % Looping through NVs to evaluate
      for nid = 1:nvnum
         n(nid) = nhs{nid}(t);
      end;
   
   % Type 3
   % - Modulating input
   case 3
      % Looping through NVs to evaluate
      for nid = 1:nvnum
         n(nid) = nhs{nid}(c(nvinps.nvminpids{nid}));
      end;
   
   % Type 4
   % - Time modulation
   % - Modulating input
   case 4
      % Looping through NVs to evaluate
      for nid = 1:nvnum
         n(nid) = nhs{nid}(t,c(nvinps.nvminpids{nid}));
      end;
   
   % Type 5
   % - Regular input
   % - Modulating input
   case 5
      % Looping through NVs to evaluate
      for nid = 1:nvnum
         n(nid) = nhs{nid}(c(nvinps.nvinpids{nid}),c(nvinps.nvminpids{nid}));
      end;
   
   % Type 6
   % - Regular input
   % - Time modulation
   case 6
      % Looping through NVs to evaluate
      for nid = 1:nvnum
         n(nid) = nhs{nid}(c(nvinps.nvinpids{nid}),t);
      end;
   
   % Type 7
   % - Regular input
   % - Time modulation
   % - Modulating input
   case 7
      % Looping through NVs to evaluate
      for nid = 1:nvnum
         n(nid) = nhs{nid}(c(nvinps.nvinpids{nid}),t,c(nvinps.nvminpids{nid}));
      end;
end;
end