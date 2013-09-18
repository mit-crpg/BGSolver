%% calcIDsForASBVs
% Right-indexing vectors for variously-sorted BVs construction function.
%
% [rASbvids,rNSbvids] = calcIDsForASBVs(srtsys,refs)
%
% This function returns two index vectors. The first indexes an
% alphabetically-sorted BV vector to order it in order of increasing BVIDs
% instead. The second reverses this operation.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       April 18, 2013
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Main function code

function [rASbvids,rNSbvids] = calcIDsForASBVs(srtsys,refs)
% This function returns the index vectors with which to index an
% alphabetically-sorted BV vector to order it in order of increasing BVIDs
% instead, and to reverse this operation. AS here means "Alphabetically-sorted",
% and NS means "Numerically-sorted" (in order of increasing BVIDs).
%
% INPUTS/OUTPUTS
% srtsys - Sorting System data structure.
% refs - Referencing arrays data structure.
% rASbvids - Index vector with which to right-index an alphabetically-sorted BV
%            vector to order it in order of increasing BVIDs instead.
% rNSbvids - Index vector with which to right-index a BV vector in order of
%            increasing BVIDs to order alphabetically instead.
% See the variable description cell in fBGSolve for more information.

% Constructing alphabetically-sorted BIDs vector of appropriate size
ASbids = calcASbids(srtsys.bnum);

% Constructing right-indexing vector for the AS BVIDs
[~,rASbvids] = sort([refs.ebvid(ASbids),refs.fbvid(ASbids)]);
rASbvids = uint16(rASbvids);

% Constructing right-indexing vector for the NS BVIDs
[~,rNSbvids] = sort(rASbvids);
rNSbvids = uint16(rNSbvids);
end

%% Alphabetically-sorted BIDs vector construction function

function ASbids = calcASbids(bnum)
% This function returns the BIDs vector of the appropriate size, with the
% elements being alphabetically-sorted.
%
% INPUTS/OUTPUTS
% bnum - Number of bonds in the system.
% ASbids - Alphabetically-sorted BIDs vector.

% Calculating max number of BID digits
dgtmax = uint8(floor(log10(single(bnum))))+1;

if dgtmax == 1
   ASbids = uint16(1:bnum);
else
   % Preallocating arrays
   ASbids = zeros(1,bnum,'uint16');
   dgtstats = false(1,bnum);
   
   % Instantiating counters
   id = uint16(1);
   ASbid = uint16(1);
   iflim = false;
   
   % Instantiating arrays
   dgtstats(1) = true;
   
   % Filling out the ASBIDs vector
   ASbids = recfillout(ASbids,dgtstats,ASbid,id,iflim,dgtmax,bnum);
end;
end

%% Recursive ASBIDs vector fillout function

function [ASbids,dgtstats,ASbid,id,iflim] = ...
   recfillout(ASbids,dgtstats,ASbid,id,iflim,dgtmax,bnum)
% This function is recursive. If the digit status vector indicates that the next
% digit to be filled is not the last digit (i.e., the resulting bid is not of
% the maximum digit length yet), a single BID entry is put into the ASbids, the
% next digit is activated, and the function recursively calls itself to continue
% the fillout. If the digit status vector indicates that the next digit to be
% filled is the last digit, then the function fills out the ASbids vector until
% reaching 9 on the last digit (or reaching the maximum BID possible), flips off
% the last digit, and returns the state. When a recursive call receives a state
% back, it increments the partial-length BID by 1 and makes the next recursive
% call, until reaching 9 on the last digit, in which case it itself returns the
% state.
%
% INPUTS/OUTPUTS
% ASbids - Alphabetically-sorted BIDs vector.
% dgtstats - Digit status vector, boolean.
% ASbid - Last alphabetically-sorted BID entered.
% id - Index in ASBIDs of the next alphabetically-sorted BID to be entered.
% iflim - Whether or not the maximum BID has been entered.
% dgtmax - Maximum number of digits in the BID.
% bnum - Number of bonds in the system.

dgtnum = uint8(nnz(dgtstats)); % Number of active digits

if dgtnum < dgtmax % Less-than-max-length ASbid
   if dgtnum == 1
      d0 = 1;
   else
      d0 = 0;
   end;
   for d = d0:9
      ASbids(id) = ASbid;
      id = id + 1;
      if (dgtnum+1 < dgtmax) || (dgtnum+1 == dgtmax && ~iflim)
         dgtstats(dgtnum+1) = true;
         ASbid = 10*ASbid;
         [ASbids,dgtstats,ASbid,id,iflim] = ...
            recfillout(ASbids,dgtstats,ASbid,id,iflim,dgtmax,bnum);
         dgtstats(dgtnum+1) = false;
      else
         ASbid = ASbid + 1;
      end;
   end;
   ASbid = uint16(floor(single(ASbid)/10));
   dgtstats(dgtnum) = false;
else % Max-length ASbid
   % DOES expect ASbid to already be evaluted to last-digit-0 form
   for d = 0:9
      ASbids(id) = ASbid;
      id = id + 1;
      ASbid = ASbid + 1;
      if ASbid > bnum
         iflim = true;
         break;
      end;
   end;
   ASbid = uint16(floor(single(ASbid)/10));
   if d < 9
      ASbid = ASbid + 1;
   end;
end;
end