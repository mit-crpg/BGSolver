%% getspatternxdot
% Function for constructing the state derivative vector's Jacobian sparsity
% pattern.
%
% srtsys = getspatternxdot(srtsys)
%
% This function constructs a double sparse matrix S of 0s and 1s. S(i,j)
% corresponds to the the partial derivative of xdot_i with respect to x_j. S is
% sized [xnum x xnum]. This matrix is the sparsity pattern of the state
% derivative vector's Jacobian.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master sparsity pattern construction function

function srtsys = getspatternxdot(srtsys)
% The purpose of this function is to construct the sparsity pattern of the state
% derivative vector's Jacobian.
%
% INPUTS/OUTPUTS
% srtsys - Sorting system data structure.

% Preallocating state variable dependency tracking arrays
xnum = double(srtsys.xnum);
xdeps = cell(1,xnum);
xdepnums = zeros(1,xnum);

% Looping through state derivative bond variables
for xid = 1:xnum
   % Retrieving state derivative BVID
   xdotID = srtsys.xdotIDs(xid);

   % Getting state derivative bond variable's layer
   xdotIDlnum = srtsys.bvlayers(xdotID);
   
   % Preallocating dependence tracking vectors
   ndstats = false(srtsys.nvnum,1);
   bdstats = false(srtsys.bvnum,1);
   bdstats(xdotID) = true;
   lastbdstats = false(srtsys.bvnum,1);
   lastbdstats(xdotID) = true;
   
   % Looping down through the state derivative BV's layers
   for lid = xdotIDlnum:(-1):1
      % Retrieving NV dependencies
      if lid > 1
         lastndstats = full(any(srtsys.bvinpnids(lastbdstats,:),1));
         ndstats(lastndstats) = true;
      end;
      
      % Retrieving BV dependencies
      if lid < xdotIDlnum
         lastbdstats = full(any(srtsys.nvdepbvids(lastndstats,:),1));
         bdstats(lastbdstats) = true;
      end;
   end;
   
   % Identifying BV and NV state dependencies
   xidsBdep = full(any(srtsys.bvinpxids(bdstats,:),1));
   xidsNdep = full(any(srtsys.nvdepxids(ndstats,:),1));
   xdeps{xid} = find(xidsBdep | xidsNdep);
   xdepnums(xid) = length(xdeps{xid});
end;

% Constructing index vectors
xRowIDs = repvectelems(1:xnum,xdepnums);
xColIDs = cell2mat(xdeps);
% Constructing the sparsity pattern matrix S
srtsys.xdotjsp = sparse(xRowIDs,xColIDs,1,xnum,xnum);
end