%% getvvids
% Variable IDs extraction function for variable vectors.
%
% [[ift,xids,bvids,nids]] = getvvids(V,vars,refs,
%     ['vars2check',vars2check,'vars2order',vars2order,'outform',outform])
%
% This function extracts variable IDs out of either a vertical vector of
% variables, or a vertical cell vector of vertical vectors of variables. The
% following variable IDs can be checked: time, state variables, bond variables
% and numeric variables. By default, time, state and bond variables are checked.
% For all variable types except for time, the order of the variables within the
% variable vector (after the unordered variables are extracted) can also be
% extracted. By default, with single vector input, the variable IDs are output
% as 2-field data structures (a field for indices vector sorted in ascending
% order and a field for variable order vector, with entries corresponding to
% indices vector). The order field is not present for unordered variables. By
% default, with cell vector of variable vectors input, the variable IDs are
% output as sparse boolean matrices (for unordered variables) and as sparse
% uint8 matrices (for ordered variables, with nonzero elements corresponding to
% position in the corresponding variable vector). Output order is the same as
% the order of the specified variables to check.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       January 15, 2013
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Master variable ID extraction function for variable vectors

function varargout = getvvids(varargin)
% The purpose of this function is to extract the corresponding variable IDs of
% variables present in the symbolic variable vector (or vectors) input to the
% function. Depending on requested outputs, not all variables are extracted.
%
% INPUTS/OUTPUTS
% V - Vertical symbolic variable vector, OR vertical cell vector of symbolic
%     variable vectors to extract variable IDs from.
% vars - Symbolic variable data structure. Contains symbolic variables that are
%        used for comparison with variables in V and extraction.
% refs - Referencing arrays data structure.
% vars2check - Variable types to extract. OPTIONAL. String, consisting of up to
%              4 characters in arbitrary order; the characters are:
%  't' - Check if time is present.
%  'x' - Check if state variables are present.
%  'b' - Check if bond variables are present.
%  'n' - Check if numeric variables are present.
%              Default value: 'txb'.
% vars2order - Whether or not to track the order of variables in the input
%              vector. OPTIONAL. Boolean vector of the same length as
%              vars2check; each element corresponds to whether or not the
%              corresponding variable type is ordered. If the vector is of wrong
%              length, it is trimmed/padded with false's. Time is always
%              overwritten to be unordered.
%              Default value: [] (all variables are unordered).
% outform - Output format descriptor. OPTIONAL. Possible values:
%  'default' - If the input is a single variable vector, 'indvect' output format
%              is used. If the input is a cell vector of variable vectors, 'mat'
%              output format is used. This is the default value of outform.
%  'inds' - Forces input tracking through either 2-field data structures (if the
%           input is a single variable vector), or through vectors of 2-field
%           data structures (if the input is a cell vector of variable vectors).
%           A single 2-field data structure has the following fields:
%   .indices - Horizontal full uint16 vector of variable IDs present in the
%              symbolic input vector. This vector is sorted in ascending order.
%   .order   - Horizontal full uint8 vector of variable positions in the
%              symbolic input vector (after unordered variables are extracted).
%              This vector's elements correspond to variable indices in
%              .indices. This field is not present for unordered variables.
%  'mat' - Forces input tracking through logical (for unordered variables) or
%          uint8 (for ordered variables) full horizontal vectors (if the input
%          is a single variable vector) or sparse matrices (if the input is a
%          cell vector of variable vectors). For ordered variables, the nonzero
%          elements in the vectors/matrices contain the corresponding variables'
%          positions in the symbolic input vector (after unordered variables are
%          extracted).
% ift - For a single variable vector input, a scalar boolean, whether or not
%       time is present in V. For a cell vector of variable vectors input, a
%       full vertical boolean vector, with elements corresponding to whether or
%       not time is present in corresponding individual variable vectors.
% xids - State variable IDs present in V, formatted according to outform.
% bvids - Bond variable IDs present in V, formatted according to outform.
% nids - Numeric variable IDs present in V, formatted according to outform.
%
% Outputs are ordered the same way vars2check is ordered.

% Input parser construction
ip = buildIP();
% Input parsing
ip.parse(varargin{:});
% Retrieving inputs
V = ip.Results.V;
vars = ip.Results.vars;
refs = ip.Results.refs;
vars2check = ip.Results.vars2check;
vars2order = ip.Results.vars2order;
outform = ip.Results.outform;
% Clearing memory
clear('ip');

% Getting variable type options
topts = getvtopts('t',vars2check,vars2order);
xopts = getvtopts('x',vars2check,vars2order);
bopts = getvtopts('b',vars2check,vars2order);
nopts = getvtopts('n',vars2check,vars2order);
vtnum = uint8(topts.ifcheck + xopts.ifcheck + bopts.ifcheck + nopts.ifcheck);

if isa(V,'sym') && isempty(V)
   ift = false;
   
   switch outform
      case {'default','inds'}
         if xopts.ifcheck
            xids.indices = zeros(1,0,'uint16');
            if xopts.iforder
               xids.order = zeros(1,0,'uint8');
            end;
         end;
         if bopts.ifcheck
            bvids.indices = zeros(1,0,'uint16');
            if bopts.iforder
               bvids.order = zeros(1,0,'uint8');
            end;
         end;
         if nopts.ifcheck
            nids.indices = zeros(1,0,'uint16');
            if nopts.iforder
               nids.order = zeros(1,0,'uint8');
            end;
         end;
      
      case 'mat'
         if xopts.ifcheck
            if xopts.iforder
               xids = zeros(1,0,'uint8');
            else
               xids = false(1,0);
            end;
         end;
         if bopts.ifcheck
            if bopts.iforder
               bvids = zeros(1,0,'uint8');
            else
               bvids = false(1,0);
            end;
         end;
         if nopts.ifcheck
            if nopts.iforder
               nids = zeros(1,0,'uint8');
            else
               nids = false(1,0);
            end;
         end;
   end;
elseif isa(V,'sym') && ~isempty(V)
   % Converting to cell vector of strings
   Vcell = cell(length(V),1);
   for vid = 1:length(Vcell)
      Vcell{vid} = char(V(vid));
   end;
   Vchar = char(Vcell);
   vchars = Vchar(:,1).';
   vchnums = Vchar(:,2:end);
   
   % Extracting unordered variables
   % Extracting time (unordered)
   if topts.ifcheck && ~topts.iforder
      tidsInVBool = (vchars == 't');
      ift = any(tidsInVBool);
      if ift
         vchars(tidsInVBool) = [];
         vchnums(tidsInVBool,:) = [];
      end;
   end;
   % Updating VID vector
   vids = uint16(str2num(vchnums)).'; %#ok<ST2NM>
   
   % Extracting state variables (unordered)
   if xopts.ifcheck && ~xopts.iforder
      qidsInVBool = (vchars == 'q');
      pidsInVBool = (vchars == 'p');
      xidsInVBool = qidsInVBool | pidsInVBool;
      if any(xidsInVBool)
         eIDsOfXIDs = vids(xidsInVBool); % EIDs of retrieved state vars
         xids.indices = sort(refs.xids(eIDsOfXIDs));
         vchars(xidsInVBool) = [];
         vids(xidsInVBool) = [];
      else
         xids.indices = zeros(1,0,'uint16');
      end;
   end;
   
   % Extracting bond variables (unordered)
   if bopts.ifcheck && ~bopts.iforder
      % Efforts
      ebvidsInVBool = (vchars == 'e');
      bIDsOfEBVIDs = vids(ebvidsInVBool); % BIDs of retrieved efforts
      % Flows
      fbvidsInVBool = (vchars == 'f');
      bIDsOfFBVIDs = vids(fbvidsInVBool); % BIDs of retrieved flows
      % Bond variables
      bvidsInVBool = ebvidsInVBool | fbvidsInVBool;
      if any(bvidsInVBool)
         bvids.indices = zeros(1,length(V),'uint16');
         bvids.indices(ebvidsInVBool) = refs.ebvid(bIDsOfEBVIDs);
         bvids.indices(fbvidsInVBool) = refs.fbvid(bIDsOfFBVIDs);
         bvids.indices = sort(bvids.indices(bvidsInVBool));
         vchars(bvidsInVBool) = [];
         vids(bvidsInVBool) = [];
      else
         bvids.indices = zeros(1,0,'uint16');
      end;
   end;
   
   % Extracting numeric variables (unordered)
   if nopts.ifcheck && ~nopts.iforder
      nidsInVBool = (vchars == 'n');
      if any(nidsInVBool)
         nids.indices = sort(vids(nidsInVBool));
         vchars(nidsInVBool) = [];
         vids(nidsInVBool) = [];
      else
         nids.indices = zeros(1,0,'uint16');
      end;
   end;
   
   % Extracting ordered variables
   % Extracting time (ordered)
   % Time is never an ordered variable
   
   % Extracting state variables (ordered)
   if xopts.ifcheck && xopts.iforder
      qidsInVBool = (vchars == 'q');
      pidsInVBool = (vchars == 'p');
      xidsInVBool = qidsInVBool | pidsInVBool;
      if any(xidsInVBool)
         eIDsOfXIDs = vids(xidsInVBool); % EIDs of retrieved state vars
         xids.indices = refs.xids(eIDsOfXIDs);
         xids.order = uint8(find(xidsInVBool));
         [xids.indices,ind] = sort(xids.indices);
         xids.order = xids.order(uint8(ind));
      else
         xids.indices = zeros(1,0,'uint16');
         xids.order = zeros(1,0,'uint8');
      end;
   end;
   
   % Extracting bond variables (ordered)
   if bopts.ifcheck && bopts.iforder
      % Efforts
      ebvidsInVBool = (vchars == 'e');
      bIDsOfEBVIDs = vids(ebvidsInVBool); % BIDs of retrieved efforts
      % Flows
      fbvidsInVBool = (vchars == 'f');
      bIDsOfFBVIDs = vids(fbvidsInVBool); % BIDs of retrieved flows
      % Bond variables
      bvidsInVBool = ebvidsInVBool | fbvidsInVBool;
      if any(bvidsInVBool)
         bvids.indices = zeros(1,length(V),'uint16');
         bvids.indices(ebvidsInVBool) = refs.ebvid(bIDsOfEBVIDs);
         bvids.indices(fbvidsInVBool) = refs.fbvid(bIDsOfFBVIDs);
         bvids.order = uint8(find(bvidsInVBool));
         [bvids.indices,ind] = sort(bvids.indices(bvidsInVBool));
         bvids.order = bvids.order(uint8(ind));
      else
         bvids.indices = zeros(1,0,'uint16');
         bvids.order = zeros(1,0,'uint8');
      end;
   end;
   
   % Extracting numeric variables (ordered)
   if nopts.ifcheck && nopts.iforder
      nidsInVBool = (vchars == 'n');
      if any(nidsInVBool)
         nids.indices = vids(nidsInVBool);
         nids.order = uint8(find(nidsInVBool));
         [nids.indices,ind] = sort(nids.indices);
         nids.order = nids.order(uint8(ind));
      else
         nids.indices = zeros(1,0,'uint16');
         nids.order = zeros(1,0,'uint8');
      end;
   end;
   
   % Converting output format
   switch outform
      case {'default','inds'}
         % Format's already correct, nothing happens
      
      case 'mat'
         if xopts.ifcheck
            if xopts.iforder
               xidsMat = zeros(1,length(vars.x),'uint8');
               xidsMat(xids.indices) = xids.order;
               xids = xidsMat;
            else
               xidsMat = false(1,length(vars.x));
               xidsMat(xids.indices) = true;
               xids = xidsMat;
            end;
         end;
         
         if bopts.ifcheck
            if bopts.iforder
               bvidsMat = zeros(1,length(vars.b),'uint8');
               bvidsMat(bvids.indices) = bvids.order;
               bvids = bvidsMat;
            else
               bvidsMat = false(1,length(vars.b));
               bvidsMat(bvids.indices) = true;
               bvids = bvidsMat;
            end;
         end;
         
         if nopts.ifcheck
            if nopts.iforder
               nidsMat = zeros(1,length(vars.n),'uint8');
               nidsMat(nids.indices) = nids.order;
               nids = nidsMat;
            else
               nidsMat = false(1,length(vars.n));
               nidsMat(nids.indices) = true;
               nids = nidsMat;
            end;
         end;
   end;
elseif iscell(V)
   % Checking output format
   if any(strcmp(outform,{'default','mat'}))
      % Counting total variables
      vnum = uint32(0);
      for veid = 1:length(V)
         vnum = vnum + length(V{veid});
      end;
      
      % Preallocating index vectors
      % Preallocating time-tracking vector
      if topts.ifcheck
         ift = false(length(V),1);
      end;
      % Preallocating state variable-tracking vectors
      if xopts.ifcheck
         xidsVids = zeros(vnum,1,'uint16');
         xidsVect = zeros(vnum,1,'uint16');
         xArrID = uint32(1);
         if xopts.iforder
            xidsInds = zeros(vnum,1,'uint8');
         end;
      end;
      % Preallocating bond variable-tracking vectors
      if bopts.ifcheck
         bvidsVids = zeros(vnum,1,'uint16');
         bvidsVect = zeros(vnum,1,'uint16');
         bvArrID = uint32(1);
         if bopts.iforder
            bvidsInds = zeros(vnum,1,'uint8');
         end;
      end;
      % Preallocating numeric variable-tracking vectors
      if nopts.ifcheck
         nidsVids = zeros(vnum,1,'uint16');
         nidsVect = zeros(vnum,1,'uint16');
         nArrID = uint32(1);
         if nopts.iforder
            nidsInds = zeros(vnum,1,'uint8');
         end;
      end;
      % Preallocating output cell array for individual subvector's tracking
      % data structures
      outstructs = cell(1,vtnum);
      
      % Looping through variable vectors
      for veid = 1:length(V)
         ve = V{veid};
         
         [outstructs{:}] = getvvids(ve,vars,refs,'vars2check',vars2check,...
            'vars2order',vars2order);
         
         % Checking if time present
         if topts.ifcheck
            ift(veid) = outstructs{topts.idInOut};
         end;
         % Checking if state variables present
         if xopts.ifcheck
            xids = outstructs{xopts.idInOut};
            xidsNum = uint32(length(xids.indices));
            xidsVids(xArrID:(xArrID+xidsNum-1)) = veid;
            xidsVect(xArrID:(xArrID+xidsNum-1)) = xids.indices;
            if xopts.iforder
               xidsInds(xArrID:(xArrID+xidsNum-1)) = xids.order;
            end;
            xArrID = xArrID + xidsNum;
         end;
         % Checking if bond variables present
         if bopts.ifcheck
            bvids = outstructs{bopts.idInOut};
            bvidsNum = uint32(length(bvids.indices));
            bvidsVids(bvArrID:(bvArrID+bvidsNum-1)) = veid;
            bvidsVect(bvArrID:(bvArrID+bvidsNum-1)) = bvids.indices;
            if bopts.iforder
               bvidsInds(bvArrID:(bvArrID+bvidsNum-1)) = bvids.order;
            end;
            bvArrID = bvArrID + bvidsNum;
         end;
         % Checking if numeric variables present
         if nopts.ifcheck
            nids = outstructs{nopts.idInOut};
            nidsNum = uint32(length(nids.indices));
            nidsVids(nArrID:(nArrID+nidsNum-1)) = veid;
            nidsVect(nArrID:(nArrID+nidsNum-1)) = nids.indices;
            if nopts.iforder
               nidsInds(nArrID:(nArrID+nidsNum-1)) = nids.order;
            end;
            nArrID = nArrID + nidsNum;
         end;
      end;
      
      % Trimming index vectors and converting to double
      if xopts.ifcheck
         xidsVids = double(nonzeros(xidsVids));
         xidsVect = double(nonzeros(xidsVect));
         if xopts.iforder
            xidsInds = double(nonzeros(xidsInds));
         end;
      end;
      if bopts.ifcheck
         bvidsVids = double(nonzeros(bvidsVids));
         bvidsVect = double(nonzeros(bvidsVect));
         if bopts.iforder
            bvidsInds = double(nonzeros(bvidsInds));
         end;
      end;
      if nopts.ifcheck
         nidsVids = double(nonzeros(nidsVids));
         nidsVect = double(nonzeros(nidsVect));
         if nopts.iforder
            nidsInds = double(nonzeros(nidsInds));
         end;
      end;
      
      % Constructing sparse matrices to return
      if xopts.ifcheck
         if xopts.iforder
            xids = sparse(xidsVids,xidsVect,xidsInds,length(V),length(vars.x));
         else
            xids = sparse(xidsVids,xidsVect,true,length(V),length(vars.x));
         end;
      end;
      if bopts.ifcheck
         if bopts.iforder
            bvids = sparse(bvidsVids,bvidsVect,bvidsInds,length(V),...
               length(vars.b));
         else
            bvids = sparse(bvidsVids,bvidsVect,true,length(V),length(vars.b));
         end;
      end;
      if nopts.ifcheck
         if nopts.iforder
            nids = sparse(nidsVids,nidsVect,nidsInds,length(V),length(vars.n));
         else
            nids = sparse(nidsVids,nidsVect,true,length(V),length(vars.n));
         end;
      end;
   elseif strcmp(outform,'inds')
      % Preallocating tracking structures
      % Preallocating time-tracking vector
      if topts.ifcheck
         ift = false(length(V),1);
      end;
      % Preallocating state variable-tracking data structure
      if xopts.ifcheck
         if xopts.iforder
            xids(1:length(V),1) = struct('indices',zeros(1,0,'uint16'),...
               'order',zeros(1,0,'uint8'));
         else
            xids(1:length(V),1) = struct('indices',zeros(1,0,'uint16'));
         end;
      end;
      % Preallocating bond variable-tracking data structure
      if bopts.ifcheck
         if bopts.iforder
            bvids(1:length(V),1) = struct('indices',zeros(1,0,'uint16'),...
               'order',zeros(1,0,'uint8'));
         else
            bvids(1:length(V),1) = struct('indices',zeros(1,0,'uint16'));
         end;
      end;
      % Preallocating numeric variable-tracking data structure
      if nopts.ifcheck
         if nopts.iforder
            nids(1:length(V),1) = struct('indices',zeros(1,0,'uint16'),...
               'order',zeros(1,0,'uint8'));
         else
            nids(1:length(V),1) = struct('indices',zeros(1,0,'uint16'));
         end;
      end;
      % Preallocating output cell array for individual subvector's tracking
      % data structures
      outstructs = cell(1,vtnum);
      
      % Looping through variable vectors
      for veid = 1:length(V)
         ve = V{veid};
         
         [outstructs{:}] = getvvids(ve,vars,refs,'vars2check',vars2check,...
            'vars2order',vars2order);
         
         % Checking if time present
         if topts.ifcheck
            ift(veid) = outstructs{topts.idInOut};
         end;
         % Checking if state variables present
         if xopts.ifcheck
            xids(veid) = outstructs{xopts.idInOut};
         end;
         % Checking if bond variables present
         if bopts.ifcheck
            bvids(veid) = outstructs{bopts.idInOut};
         end;
         % Checking if numeric variables present
         if nopts.ifcheck
            nids(veid) = outstructs{nopts.idInOut};
         end;
      end;
   end;
end;

% Constructing output array
% Preallocating output array
varargout = cell(1,topts.ifcheck+xopts.ifcheck+bopts.ifcheck+nopts.ifcheck);
% Recording if time present
if topts.ifcheck
   varargout{topts.idInOut} = ift;
end;
% Recording state variable IDs
if xopts.ifcheck
   varargout{xopts.idInOut} = xids;
end;
% Recording bond variable IDs
if bopts.ifcheck
   varargout{bopts.idInOut} = bvids;
end;
% Recording numeric variable IDs
if nopts.ifcheck
   varargout{nopts.idInOut} = nids;
end;
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for
% the getvvids function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the getvvids function.

fVcheck = @(V) isvector(V) && (iscell(V) || isa(V,'sym'));
fVars2Check = @(vars2check) all(ismember(vars2check,'txbn'));
fVars2Order = @(vars2order) isempty(vars2order) || islogical(vars2order);
fOutFormCheck = @(retform) any(strcmp(retform,{'default','mat','inds'}));

ip = inputParser;
ip.addRequired('V',fVcheck);
ip.addRequired('vars',@isstruct);
ip.addRequired('refs',@isstruct);
ip.addParamValue('vars2check','txb',fVars2Check);
ip.addParamValue('vars2order',false(0),fVars2Order);
ip.addParamValue('outform','default',fOutFormCheck);
end