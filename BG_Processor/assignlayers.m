%% assignlayers
% Numeric layer assignment function for bond and numeric variables.
%
% srtsys = assignlayers(srtsys,runopts,refs)
%
% This function assigns numeric layers to all bond and numeric variables in the
% BGS. It takes a sorting system data structure that contains the bond variable
% expressions vector and the bond and numeric variable input tracking arrays.
% The function appends numeric layer tracking arrays to the data structure and
% returns it. Numeric layers work as follows: bond variables on numeric layer n
% have to be evaluated before numeric variables on numeric layer n can be
% evaluated. Bond and numeric variables on numeric layer n have to be evaluated
% before bond and numeric variables on numeric layers >n can be evaluated.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Main numeric layer assignment function code

function srtsys = assignlayers(srtsys)
% The purpose of this function is to assign numeric layers to all bond and
% numeric variables in the BGS. It takes a sorting system data structure that
% contains the bond variable expressions vector and the bond and numeric
% variable input tracking arrays. The function appends numeric layer tracking
% arrays to the data structure and returns it.
%
% INPUTS/OUTPUTS
% srtsys - Sorting system data structure, with bond variable expressions vector
%          and the bond and numeric variable input tracking arrays already
%          augmented to it. The function assumes the input tracking arrays to be
%          in sparse boolean/double matrices format.

% Preallocating arrays
srtsys.bvlayers = zeros(1,srtsys.bvnum,'uint8');
srtsys.nvlayers = zeros(1,srtsys.nvnum,'uint8');
bvlstat = false(srtsys.bvnum,1);
nvlstat = false(srtsys.nvnum,1);

% Instantiating counters
srtsys.lnum = uint8(0); % Number of numeric layers

% Looping through layers, until all variables have a layer assigned
while ~all([bvlstat;nvlstat])
   % Incrementing layer counter
   srtsys.lnum = srtsys.lnum + 1;
   
   % Identifying bond variables on this layer
   bvids = false(srtsys.bvnum,1);
   bvids(~bvlstat) = full(~any(srtsys.bvinpnids(~bvlstat,~nvlstat),2));
   srtsys.bvlayers(bvids) = srtsys.lnum;
   bvlstat(bvids) = true;
   
   % Identifying numeric variables on this layer
   nvids = false(srtsys.nvnum,1);
   nvids(~nvlstat) = full(~any(srtsys.nvdepbvids(~nvlstat,~bvlstat),2));
   srtsys.nvlayers(nvids) = srtsys.lnum;
   nvlstat(nvids) = true;
end;

% Constructing BV numeric layer tracking arrays
srtsys.layerbvs = cell(srtsys.lnum,1);
for lid = 1:srtsys.lnum
   bvsl = (srtsys.bvlayers == lid);
   srtsys.layerbvs{lid} = uint16(find(bvsl));
end;

% Constructing NV numeric layer tracking arrays
srtsys.layernvs = cell(srtsys.lnum-1,1);
for lid = 1:srtsys.lnum-1
   srtsys.layernvs{lid} = uint16(find(srtsys.nvlayers == lid));
end;
end