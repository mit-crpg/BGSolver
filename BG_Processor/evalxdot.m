%% evalxdot
% State derivative vector evaluation function MATLAB .m file.
%
% xdot = evalxdot(t,x,xdotsys)
%
% This function evaluates the state derivative vector at a given time, with a
% given (possibly vectorized) x vector, based on a given state derivative vector
% system data structure. It is used to evaluate the state derivative vector if a
% compiled .mex version is not available for the platform MATLAB is currently
% running on.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% State derivative vector evaluation function

function xdot = evalxdot(t,x,xdotsys)
% The purpose of this function is to evaluate the state derivative vector at a
% given time, with a given (possibly vectorized) x vector, based on a given
% state derivative vector system data structure.
%
% INPUTS/OUTPUTS
% t - Time to input to xdot. double scalar.
% x - State vector to input to xdot, possibly vectorized such that x(:,m) is the
%     m'th vertical state vector to input to xdot. Vertical double vector OR a
%     matrix of vertical double vectors.
% xdotsys - State derivative vector system data structure, created by
%           createxdotsys function.
% xdot - State derivative vector evaluated at time t and with state vector x. If
%        x was vectorized, xdot is a matrix, such that xdot(:,m) is the m'th
%        vertical state derivative vector at time t and state x(:,m).

% Sizing the problem
nvect = uint32(size(x,2));

% Preallocating arrays
brfx = zeros(xdotsys.bvnumrfx,nvect);
nrfx = zeros(xdotsys.nvnumrfx,nvect);

% Going through every layer and evaluating bond and numeric variables
for lid = 1:xdotsys.lnumrfx
   % Evaluating bond variables on the layer
   if ~isempty(xdotsys.lbvidsrfx{lid})
      bvinvars = cell(1,xdotsys.bvInpSiz(lid));
      if xdotsys.bvifit(lid)
         bvinvars{xdotsys.bvtPosInInp(lid)} = t;
      end;
      if xdotsys.bvifix(lid)
         bvinvars{xdotsys.bvxPosInInp(lid)} = x(xdotsys.lbvinpxids{lid},:);
      end;
      if xdotsys.bvifin(lid)
         bvinvars{xdotsys.bvnPosInInp(lid)} = ...
            nrfx(xdotsys.lbvinpnidsrfx{lid},:);
      end;
      if lid == 1 && xdotsys.bv1vectstat == 2
         scallbvvect = xdotsys.b{lid}(bvinvars{:});
         brfx(xdotsys.lbvidsrfx{lid},:) = repmat(scallbvvect,1,nvect);
      else
         brfx(xdotsys.lbvidsrfx{lid},:) = xdotsys.b{lid}(bvinvars{:});
      end;
   end;
   
   % Evaluating numeric variables on the layer
   if lid < xdotsys.lnumrfx
      % Constructing concatenated input vector
      c = [x;brfx];
      % Looping through NV types
      for nvtyp = 1:7
         tnidsrfx = xdotsys.lnvidsrfxByType{nvtyp,lid};
         tnhs = xdotsys.nhs(tnidsrfx);
         tnvinps = xdotsys.lnvinpsByType{nvtyp,lid};
         nrfx(tnidsrfx,:) = evalnv(t,c,nvtyp,tnhs,tnvinps);
      end;
   end;
end;

% Constructing the state derivative vector
xdot = brfx(xdotsys.xdotIDsrfx,:);
end

%% Numeric variable evaluation function

function n = evalnv(t,c,nvtyp,nhs,nvinps)
% The purpose of this function is to evaluate a vector of numeric variables of a
% given type.
%
% INPUTS/OUTPUTS
% t - Time to evaluate the numeric variables at.
% c - Concatenated input vector/matrix.
% nidsrfx - NIDs RFX of the numeric variables to evaluate.
% nvtype - Numeric variable type of the numeric variables to evaluate. See
%          createxdotsys variable documentation cell for details of NV types.
% nhs - Vector of numeric function handles to evaluate.
% nvinps - NV input tracking data structure.
% n - Evaluated numeric variable vector/matrix.

% Sizing the problem
nvnum = uint16(length(nhs));

% Prealloacting n
if nvtyp == 2
   n = zeros(nvnum,1);
elseif any(nvtyp == [1,3:7])
   n = zeros(nvnum,size(c,2));
end;

switch nvtyp
   % Type 1
   % - Regular input
   case 1
      % Looping through NVs to evaluate
      for nid = 1:nvnum
         n(nid,:) = nhs{nid}(c(nvinps.nvinpids{nid},:));
      end;
   
   % Type 2
   % - Time modulation
   case 2
      % Looping through NVs to evaluate
      for nid = 1:nvnum
         n(nid) = nhs{nid}(t);
      end;
      % Vectorizing the output
      n = repmat(n,1,size(c,2));
   
   % Type 3
   % - Modulating input
   case 3
      % Looping through NVs to evaluate
      for nid = 1:nvnum
         n(nid,:) = nhs{nid}(c(nvinps.nvminpids{nid},:));
      end;
   
   % Type 4
   % - Time modulation
   % - Modulating input
   case 4
      % Looping through NVs to evaluate
      for nid = 1:nvnum
         n(nid,:) = nhs{nid}(t,c(nvinps.nvminpids{nid},:));
      end;
   
   % Type 5
   % - Regular input
   % - Modulating input
   case 5
      % Looping through NVs to evaluate
      for nid = 1:nvnum
         n(nid,:) = nhs{nid}(c(nvinps.nvinpids{nid},:),...
            c(nvinps.nvminpids{nid},:));
      end;
   
   % Type 6
   % - Regular input
   % - Time modulation
   case 6
      % Looping through NVs to evaluate
      for nid = 1:nvnum
         n(nid,:) = nhs{nid}(c(nvinps.nvinpids{nid},:),t);
      end;
   
   % Type 7
   % - Regular input
   % - Time modulation
   % - Modulating input
   case 7
      % Looping through NVs to evaluate
      for nid = 1:nvnum
         n(nid,:) = nhs{nid}(c(nvinps.nvinpids{nid},:),t,...
               c(nvinps.nvminpids{nid},:));
      end;
end;
end