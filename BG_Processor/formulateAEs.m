%% formulateAEs
% Algebraic equation formulation function.
%
% srtsys = formulateAEs(bgs,refs,vars,ebcm)
%
% This function formulates the algebraic equations for a given bond graph
% system. It accepts the bond graph system data structure, the referencing
% arrays data structure and the symbolic variables data structure. The algebraic
% equations are returned in a Sorting System data structure, along with
% supplementary indexing information.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Variable descriptions
%
% INPUTS
% - All inputs are required.
%
% bgs - Bond Graph System data structure. Describes sections 7-9 in the BGSD
%       file. It is assumed that the data structure has already been completed
%       by completeBGSDdata, with BGSolve operation mode. Therefore, assuming EL
%       storage type a.
% refs - Referencing arrays data structure. Contains index arrays that are used
%        to reference quantities by EID (or BID), and not by quantity ID.
%        refs.Quantity(EID) is the quantity ID for quantity Quantity of element
%        EID.
% vars - Symbolic variable data structure. Contains symbolic variables that are
%        used for sorting and formulation of algebraic equations.
% ebcm - Element-indexed Bond Connectivity Map data structure.
%
% OUTPUTS
% srtsys - Sorting System data structure.
%
% See the variable description cell in fBGSolve function for more information.

%% Master algebraic equation formulation function

function srtsys = formulateAEs(varargin)
% The purpose of this function is to formulate the algebraic equations for a
% given bond graph system and return them in a sorting system data structure.
%
% INPUTS/OUTSPUTS
% bgs - Bond Graph System data structure. EL storage type a.
% refs - Referencing arrays data structure.
% vars - Variables data structure.
% srtsys - Sorting System data structure.
% See the variable description cell for more information.

% Constructing input parser
ip = buildIP();
% Parsing input
ip.parse(varargin{:});
% Retrieving inputs
bgs = ip.Results.bgs;
refs = ip.Results.refs;
vars = ip.Results.vars;
ebcm = ip.Results.ebcm;
% Clearing memory
clear('ip');

% Sizing arrays
srtsys.bvnum = 2*bgs.bnum;
srtsys.bnum = bgs.bnum;
srtsys.xnum = bgs.xnum;
srtsys.nvnum = bgs.nvnum;

% Preallocating and copying arrays
srtsys.aes = cell(srtsys.bvnum,1);
srtsys.xdotIDs = zeros(1,bgs.xnum,'uint16');
srtsys.xphxs = zeros(bgs.xnum,1,'uint8');
srtsys.xnns = false(bgs.xnum,1);
srtsys.xats = zeros(bgs.xnum,1);
srtsys.nhs = cell(bgs.nvnum,1);
srtsys.nvinpxids(1:bgs.nvnum,1) = struct('indices',zeros(1,0,'uint16'),...
   'order',zeros(1,0,'uint8')); % XID input tracking; could be a uint16 vector,
                                % but implemented this way for extensibility
srtsys.nvinpbvids(1:bgs.nvnum,1) = struct('indices',zeros(1,0,'uint16'),...
   'order',zeros(1,0,'uint8')); % BVID input tracking
srtsys.nvmitids = false(bgs.nvnum,1); % Time modulation tracking
srtsys.nvminpxids(1:bgs.nvnum,1) = struct('indices',zeros(1,0,'uint16'),...
   'order',zeros(1,0,'uint8')); % Modulating XID input tracking
srtsys.nvminpbvids(1:bgs.nvnum,1) = struct('indices',zeros(1,0,'uint16'),...
   'order',zeros(1,0,'uint8')); % Modulating BVID input tracking
srtsys.x0s = bgs.x0s;

% Instantiating counters
aeid = uint16(1);
nid = uint16(1);

% Formulating algebraic equations for every element
for eid = 1:bgs.enum
   % Extracting element type
   etype = bgs.etypes{eid};
   
   % Formulating element-specific algebraic equation(s)
   switch etype
      % Source of effort
      % Assuming:
      % - exprype = CC
      case 'SE'
         bid = ebcm.conn{eid};
         exprtype = bgs.exprtypes{refs.exprs(eid)};
         expr = bgs.exprs{refs.exprs(eid)};
         
         % Formulating the constituent equation
         % e = const
         if strcmp(exprtype,'CC')
            srtsys.aes{aeid} = vars.e(bid) - expr;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 1;
      
      % Source of flow
      % Assuming:
      % - exprtype = CC
      case 'SF'
         bid = ebcm.conn{eid};
         exprtype = bgs.exprtypes{refs.exprs(eid)};
         expr = bgs.exprs{refs.exprs(eid)};
         
         % Formulating the constituent equation
         % f = const
         if strcmp(exprtype,'CC')
            srtsys.aes{aeid} = vars.f(bid) - expr;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 1;
      
      % Inertial element
      % Assuming:
      % - exprtype = CC, CE or NE
      case 'I'
         bid = ebcm.conn{eid};
         xid = refs.xids(eid);
         exprid = refs.exprs(eid);
         phxid = refs.phxs(eid);
         nnid = refs.nns(eid);
         atid = refs.ats(eid);
         srtsys.xdotIDs(xid) = refs.ebvid(bid);
         srtsys.xphxs(xid) = bgs.phxs(phxid);
         srtsys.xnns(xid) = bgs.nns(nnid);
         srtsys.xats(xid) = bgs.ats(atid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         
         % Formulating the constituent equation
         % f = I(p)
         if strcmp(exprtype,'CC')
            srtsys.aes{aeid} = vars.f(bid) - expr * vars.x(xid);
         elseif strcmp(exprtype,'CE')
            exprSub = subs(expr,vars.ui.p,vars.x(xid));
            srtsys.aes{aeid} = vars.f(bid) - exprSub;
         elseif strcmp(exprtype,'NE')
            % Recording numeric indices
            srtsys.nvinpxids(nid) = getinpids(xid);
            % Constructing equation
            srtsys.nhs{nid} = expr;
            srtsys.aes{aeid} = vars.f(bid) - vars.n(nid);
            % Incrementing indices
            nid = nid + 1;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 1;
      
      % Capacitive element
      % Assuming:
      % - exprtype = CC, CE or NE
      case 'C'
         bid = ebcm.conn{eid};
         xid = refs.xids(eid);
         exprid = refs.exprs(eid);
         phxid = refs.phxs(eid);
         nnid = refs.nns(eid);
         atid = refs.ats(eid);
         srtsys.xdotIDs(xid) = refs.fbvid(bid);
         srtsys.xphxs(xid) = bgs.phxs(phxid);
         srtsys.xnns(xid) = bgs.nns(nnid);
         srtsys.xats(xid) = bgs.ats(atid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         
         % Formulating the constituent equation
         % e = C(p)
         if strcmp(exprtype,'CC')
            srtsys.aes{aeid} = vars.e(bid) - expr * vars.x(xid);
         elseif strcmp(exprtype,'CE')
            exprSub = subs(expr,vars.ui.q,vars.x(xid));
            srtsys.aes{aeid} = vars.e(bid) - exprSub;
         elseif strcmp(exprtype,'NE')
            % Recording numeric indices
            srtsys.nvinpxids(nid) = getinpids(xid);
            % Constructing equation
            srtsys.nhs{nid} = expr;
            srtsys.aes{aeid} = vars.e(bid) - vars.n(nid);
            % Incrementing indices
            nid = nid + 1;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 1;
      
      % Resistive element
      % Assuming:
      % - exprtype = CC, CE or NE
      case 'R'
         bid = ebcm.conn{eid};
         exprid = refs.exprs(eid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         
         % Formulating the constituent equation
         % e = R(f)
         %  or
         % f = G(e)
         if strcmp(exprtype,'CC')
            srtsys.aes{aeid} = vars.e(bid) - expr * vars.f(bid);
         elseif strcmp(exprtype,'CE')
            % Retrieving causality
            cambid = refs.causes(eid);
            caus = bgs.causes{cambid};
            % Checking causality
            if caus == 1 % Effort is input, flow is output
               ivarUI = vars.ui.e;
               ivar = vars.e(bid);
               ovar = vars.f(bid);
            elseif caus == 2 % Flow is input, effort is output
               ivarUI = vars.ui.f;
               ivar = vars.f(bid);
               ovar = vars.e(bid);
            else
               error('EID %u, ET %s, EXPRT %s, unexpected EC!!!',eid,etype,...
                  exprtype);
            end;
            % Constructing equation
            srtsys.aes{aeid} = ovar - subs(expr,ivarUI,ivar);
         elseif strcmp(exprtype,'NE')
            % Retrieving causality
            cambid = refs.causes(eid);
            caus = bgs.causes{cambid};
            % Checking causality
            if caus == 1 % Effort is input, flow is output
               ovar = vars.f(bid);
               % Recording numeric indices
               srtsys.nvinpbvids(nid) = getinpids(refs.ebvid(bid));
            elseif caus == 2 % Flow is input, effort is output
               ovar = vars.e(bid);
               % Recording numeric indices
               srtsys.nvinpbvids(nid) = getinpids(refs.fbvid(bid));
            else
               error('EID %u, ET %s, EXPRT %s, unexpected EC!!!',eid,etype,...
                  exprtype);
            end;
            % Constructing equation
            srtsys.nhs{nid} = expr;
            srtsys.aes{aeid} = ovar - vars.n(nid);
            % Incrementing indices
            nid = nid + 1;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 1;
      
      % Transformer element
      % Assuming:
      % - exprtype = CC
      % - Exactly two bonds connected to the element
      % - One bond points to the element
      % - One bond points from the element
      case 'TF'
         bidFrom = ebcm.from{eid};
         bidTo = ebcm.to{eid};
         exprid = refs.exprs(eid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         
         % Formulating the constituent equations
         % e_bidTo   = m * e_bidFrom
         % f_bidFrom = m * f_bidTo
         if strcmp(exprtype,'CC')
            srtsys.aes{aeid} = vars.e(bidTo) - expr * vars.e(bidFrom);
            srtsys.aes{aeid+1} = vars.f(bidFrom) - expr * vars.f(bidTo);
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 2;
      
      % Gyrator element
      % Assuming:
      % - exprtype = CC
      % - Exactly two bonds connected to the element
      % - One bond points to the element
      % - One bond points from the element
      case 'GY'
         bidFrom = ebcm.from{eid};
         bidTo = ebcm.to{eid};
         exprid = refs.exprs(eid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         
         % Formulating the constituent equations
         % e_bidTo   = r * f_bidFrom
         % e_bidFrom = r * f_bidTo
         if strcmp(exprtype,'CC')
            srtsys.aes{aeid} = vars.e(bidTo) - expr * vars.f(bidFrom);
            srtsys.aes{aeid+1} = vars.e(bidFrom) - expr * vars.f(bidTo);
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 2;
      
      % 1-junction element
      % Assuming:
      % - Two or more bonds are connected to the junction
      %  OR
      % - A single bond is connected to the junction, with flow causality
      %   towards the junction
      case '1'
         bidsConn = ebcm.conn{eid};
         bidsFrom = ebcm.from{eid};
         bidsTo = ebcm.to{eid};
         bcnum = uint16(length(bidsConn));
         
         % Formulating the constituent equations
         if bcnum == 1
            % Effort equation
            srtsys.aes{aeid} = vars.e(bidsConn);
            aeid = aeid + 1;
         else
            % Flow equations
            % f_i = f_i+1 = ... = f_i+k
            ivect = ones(bcnum-1,1) * vars.f(bidsConn(1));
            ovect = vars.f(bidsConn(2:end));
            srtsys.aes(aeid:(aeid+bcnum-2)) = num2cell(ivect - ovect);
            aeid = aeid + bcnum - 1;
            % Effort equations
            % sum(e_to) - sum(e_from) = 0
            srtsys.aes{aeid} = sum(vars.e(bidsTo)) - sum(vars.e(bidsFrom));
            aeid = aeid + 1;
         end;
      
      % 0-junction element
      % Assuming:
      % - Two or more bonds are connected to the junction
      %  OR
      % - A single bond is connected to the junction, with effort causality
      %   towards the junction
      case '0'
         bidsConn = ebcm.conn{eid};
         bidsFrom = ebcm.from{eid};
         bidsTo = ebcm.to{eid};
         bcnum = uint16(length(bidsConn));
         
         % Formulating the constituent equations
         if bcnum == 1
            % Flow equation
            srtsys.aes{aeid} = vars.f(bidsConn);
            aeid = aeid + 1;
         else
            % Effort equations
            % e_i = e_i+1 = ... = e_i+k
            ivect = ones(bcnum-1,1) * vars.e(bidsConn(1));
            ovect = vars.e(bidsConn(2:end));
            srtsys.aes(aeid:(aeid+bcnum-2)) = num2cell(ivect - ovect);
            aeid = aeid + bcnum - 1;
            % Flow equations
            % sum(f_to) - sum(f_from) = 0
            srtsys.aes{aeid} = sum(vars.f(bidsTo)) - sum(vars.f(bidsFrom));
            aeid = aeid + 1;
         end;
      
      % Modulated source of effort
      % Assuming:
      % - exprtype = CMC or NMC
      case 'MSE'
         bid = ebcm.conn{eid};
         exprid = refs.exprs(eid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         
         % Formulating the constituent equation
         % e = func_m([t],[m])
         if strcmp(exprtype,'CMC')
            srtsys.aes{aeid} = vars.e(bid) - expr;
         elseif strcmp(exprtype,'NMC')
            % Retrieving modulating variable list
            modid = refs.mvars(eid);
            mvars = bgs.mvars{modid};
            % Identifying modulating variables
            [srtsys.nvmitids(nid),srtsys.nvminpxids(nid),...
               srtsys.nvminpbvids(nid)] = getvvids(mvars,vars,refs,...
               'vars2check','txb','vars2order',[false,true,true]);
            % Constructing equation
            srtsys.nhs{nid} = expr;
            srtsys.aes{aeid} = vars.e(bid) - vars.n(nid);
            % Incrementing indices
            nid = nid + 1;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 1;
      
      % Modulated source of flow
      % Assuming:
      % - exprtype = CMC or NMC
      case 'MSF'
         bid = ebcm.conn{eid};
         exprid = refs.exprs(eid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         
         % Formulating the constituent equation
         % f = func_m([t],[m])
         if strcmp(exprtype,'CMC')
            srtsys.aes{aeid} = vars.f(bid) - expr;
         elseif strcmp(exprtype,'NMC')
            % Retrieving modulating variable list
            modid = refs.mvars(eid);
            mvars = bgs.mvars{modid};
            % Identifying modulating variables
            [srtsys.nvmitids(nid),srtsys.nvminpxids(nid),...
               srtsys.nvminpbvids(nid)] = getvvids(mvars,vars,refs,...
               'vars2check','txb','vars2order',[false,true,true]);
            % Constructing equation
            srtsys.nhs{nid} = expr;
            srtsys.aes{aeid} = vars.f(bid) - vars.n(nid);
            % Incrementing indices
            nid = nid + 1;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 1;
      
      % Modulated inertial element
      % Assuming:
      % - exprtype = CMC, NMC, CME or NME
      case 'MI'
         bid = ebcm.conn{eid};
         xid = refs.xids(eid);
         exprid = refs.exprs(eid);
         nnid = refs.nns(eid);
         atid = refs.ats(eid);
         srtsys.xdotIDs(xid) = refs.ebvid(bid);
         srtsys.xnns(xid) = bgs.nns(nnid);
         srtsys.xats(xid) = bgs.ats(atid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         
         % Formulating the constituent equation
         % f = I_m(p,[t],[m])
         if strcmp(exprtype,'CMC')
            srtsys.aes{aeid} = vars.f(bid) - expr * vars.x(xid);
         elseif strcmp(exprtype,'NMC')
            % Retrieving modulating variable list
            modid = refs.mvars(eid);
            mvars = bgs.mvars{modid};
            % Identifying modulating variables
            [srtsys.nvmitids(nid),srtsys.nvminpxids(nid),...
               srtsys.nvminpbvids(nid)] = getvvids(mvars,vars,refs,...
               'vars2check','txb','vars2order',[false,true,true]);
            % Constructing equation
            srtsys.nhs{nid} = expr;
            srtsys.aes{aeid} = vars.f(bid) - vars.n(nid) * vars.x(xid);
            % Incrementing indices
            nid = nid + 1;
         elseif strcmp(exprtype,'CME')
            exprSub = subs(expr,vars.ui.p,vars.x(xid));
            srtsys.aes{aeid} = vars.f(bid) - exprSub;
         elseif strcmp(exprtype,'NME')
            % Recording numeric indices
            srtsys.nvinpxids(nid) = getinpids(xid);
            % Retrieving modulating variable list
            modid = refs.mvars(eid);
            mvars = bgs.mvars{modid};
            % Identifying modulating variables
            [srtsys.nvmitids(nid),srtsys.nvminpxids(nid),...
               srtsys.nvminpbvids(nid)] = getvvids(mvars,vars,refs,...
               'vars2check','txb','vars2order',[false,true,true]);
            % Constructing equation
            srtsys.nhs{nid} = expr;
            srtsys.aes{aeid} = vars.f(bid) - vars.n(nid);
            % Incrementing indices
            nid = nid + 1;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 1;
      
      % Modulated capacitive element
      % Assuming:
      % - exprtype = CMC, NMC, CME or NME
      case 'MC'
         bid = ebcm.conn{eid};
         xid = refs.xids(eid);
         exprid = refs.exprs(eid);
         nnid = refs.nns(eid);
         atid = refs.ats(eid);
         srtsys.xdotIDs(xid) = refs.fbvid(bid);
         srtsys.xnns(xid) = bgs.nns(nnid);
         srtsys.xats(xid) = bgs.ats(atid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         
         % Formulating the constituent equation
         % e = C_m(q,[t],[m])
         if strcmp(exprtype,'CMC')
            srtsys.aes{aeid} = vars.e(bid) - expr * vars.x(xid);
         elseif strcmp(exprtype,'NMC')
            % Retrieving modulating variable list
            modid = refs.mvars(eid);
            mvars = bgs.mvars{modid};
            % Identifying modulating variables
            [srtsys.nvmitids(nid),srtsys.nvminpxids(nid),...
               srtsys.nvminpbvids(nid)] = getvvids(mvars,vars,refs,...
               'vars2check','txb','vars2order',[false,true,true]);
            % Constructing equation
            srtsys.nhs{nid} = expr;
            srtsys.aes{aeid} = vars.e(bid) - vars.n(nid) * vars.x(xid);
            % Incrementing indices
            nid = nid + 1;
         elseif strcmp(exprtype,'CME')
            exprSub = subs(expr,vars.ui.q,vars.x(xid));
            srtsys.aes{aeid} = vars.e(bid) - exprSub;
         elseif strcmp(exprtype,'NME')
            % Recording numeric indices
            srtsys.nvinpxids(nid) = getinpids(xid);
            % Retrieving modulating variable list
            modid = refs.mvars(eid);
            mvars = bgs.mvars{modid};
            % Identifying modulating variables
            [srtsys.nvmitids(nid),srtsys.nvminpxids(nid),...
               srtsys.nvminpbvids(nid)] = getvvids(mvars,vars,refs,...
               'vars2check','txb','vars2order',[false,true,true]);
            % Constructing equation
            srtsys.nhs{nid} = expr;
            srtsys.aes{aeid} = vars.f(bid) - vars.n(nid);
            % Incrementing indices
            nid = nid + 1;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 1;
      
      % Modulated resistive element
      % Assuming:
      % - exprType = CMC, NMC, CME or NME
      case 'MR'
         bid = ebcm.conn{eid};
         exprid = refs.exprs(eid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         
         % Formulating the constituent equation
         % e = R(f,[t],[m])
         %  or
         % f = G(e,[t],[m])
         if strcmp(exprtype,'CMC')
            srtsys.aes{aeid} = vars.e(bid) - expr * vars.f(bid);
         elseif strcmp(exprtype,'NMC')
            % Retrieving modulating variable list
            modid = refs.mvars(eid);
            mvars = bgs.mvars{modid};
            % Identifying modulating variables
            [srtsys.nvmitids(nid),srtsys.nvminpxids(nid),...
               srtsys.nvminpbvids(nid)] = getvvids(mvars,vars,refs,...
               'vars2check','txb','vars2order',[false,true,true]);
            % Constructing equation
            srtsys.nhs{nid} = expr;
            srtsys.aes{aeid} = vars.e(bid) - vars.n(nid) * vars.f(bid);
            % Incrementing indices
            nid = nid + 1;
         elseif strcmp(exprtype,'CME')
            % Retrieving causality
            cambid = refs.causes(eid);
            caus = bgs.causes{cambid};
            % Checking causality
            if caus == 1 % Effort is input, flow is output
               ivarUI = vars.ui.e;
               ivar = vars.e(bid);
               ovar = vars.f(bid);
            elseif caus == 2 % Flow is input, effort is output
               ivarUI = vars.ui.f;
               ivar = vars.f(bid);
               ovar = vars.e(bid);
            else
               error('EID %u, ET %s, EXPRT %s, unexpected EC!!!',eid,etype,...
                  exprtype);
            end;
            % Constructing equation
            srtsys.aes{aeid} = ovar - subs(expr,ivarUI,ivar);
         elseif strcmp(exprtype,'NME')
            % Retrieving modulating variable list
            modid = refs.mvars(eid);
            mvars = bgs.mvars{modid};
            % Identifying modulating variables
            [srtsys.nvmitids(nid),srtsys.nvminpxids(nid),...
               srtsys.nvminpbvids(nid)] = getvvids(mvars,vars,refs,...
               'vars2check','txb','vars2order',[false,true,true]);
            % Retrieving causality
            cambid = refs.causes(eid);
            caus = bgs.causes{cambid};
            % Checking causality
            if caus == 1 % Effort is input, flow is output
               ovar = vars.f(bid);
               % Recording numeric indices
               srtsys.nvinpbvids(nid) = getinpids(refs.ebvid(bid));
            elseif caus == 2 % Flow is input, effort is output
               ovar = vars.e(bid);
               % Recording numeric indices
               srtsys.nvinpbvids(nid) = getinpids(refs.fbvid(bid));
            else
               error('EID %u, ET %s, EXPRT %s, unexpected EC!!!',eid,etype,...
                  exprtype);
            end;
            % Constructing equation
            srtsys.nhs{nid} = expr;
            srtsys.aes{aeid} = ovar - vars.n(nid);
            % Incrementing indices
            nid = nid + 1;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 1;
      
      % Modulated transformer element
      % Assuming:
      % - exprtype = CMC or NMC
      % - Exactly two bonds connected to the element
      % - One bond points to the element
      % - One bond points from the element
      case 'MTF'
         bidFrom = ebcm.from{eid};
         bidTo = ebcm.to{eid};
         exprid = refs.exprs(eid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         
         % Formulating the constituent equations
         % e_bidTo   = m([t],[m]) * e_bidFrom
         % f_bidFrom = m([t],[m]) * f_bidTo
         if strcmp(exprtype,'CMC')
            srtsys.aes{aeid} = vars.e(bidTo) - expr * vars.e(bidFrom);
            srtsys.aes{aeid+1} = vars.f(bidFrom) - expr * vars.f(bidTo);
         elseif strcmp(exprtype,'NMC')
            % Retrieving modulating variable list
            modid = refs.mvars(eid);
            mvars = bgs.mvars{modid};
            % Identifying modulating variables
            [srtsys.nvmitids(nid),srtsys.nvminpxids(nid),...
               srtsys.nvminpbvids(nid)] = getvvids(mvars,vars,refs,...
               'vars2check','txb','vars2order',[false,true,true]);
            % Constructing equations
            srtsys.nhs{nid} = expr;
            srtsys.aes{aeid} = vars.e(bidTo) - vars.n(nid) * vars.e(bidFrom);
            srtsys.aes{aeid+1} = vars.f(bidFrom) - vars.n(nid) * vars.f(bidTo);
            % Incrementing indices
            nid = nid + 1;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 2;
      
      % Modulating gyrator element
      % Assuming:
      % - exprtype = CMC or NMC
      % - Exactly two bonds connected to the element
      % - One bond points to the element
      % - One bond points from the element
      case 'MGY'
         bidFrom = ebcm.from{eid};
         bidTo = ebcm.to{eid};
         exprid = refs.exprs(eid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         
         % Formulating the constituent equations
         % e_bidTo   = r([t],[m]) * f_bidFrom
         % e_bidFrom = r([t],[m]) * f_bidTo
         if strcmp(exprtype,'CMC')
            srtsys.aes{aeid} = vars.e(bidTo) - expr * vars.f(bidFrom);
            srtsys.aes{aeid+1} = vars.e(bidFrom) - expr * vars.f(bidTo);
         elseif strcmp(exprtype,'NMC')
            % Retrieving modulating variable list
            modid = refs.mvars(eid);
            mvars = bgs.mvars{modid};
            % Identifying modulating variables
            [srtsys.nvmitids(nid),srtsys.nvminpxids(nid),...
               srtsys.nvminpbvids(nid)] = getvvids(mvars,vars,refs,...
               'vars2check','txb','vars2order',[false,true,true]);
            % Constructing equations
            srtsys.nhs{nid} = expr;
            srtsys.aes{aeid} = vars.e(bidTo) - vars.n(nid) * vars.f(bidFrom);
            srtsys.aes{aeid+1} = vars.e(bidFrom) - vars.n(nid) * vars.f(bidTo);
            % Incrementing indices
            nid = nid + 1;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 2;
         
      % Nonconservative 2-port coupler element
      % Assuming:
      % - exprtype = CC, CE or NE
      % - Exactly two bonds are connected to the element
      % - One bond points to the element
      % - One bond points from the element
      case 'R2'
         bidFrom = ebcm.from{eid};
         bidTo = ebcm.to{eid};
         exprid = refs.exprs(eid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         cambid = refs.causes(eid);
         caus = bgs.causes{cambid};
         
         % Checking causality
         if all(caus == [1,1]) % Inputs: e_i, e_o | Outputs: f_i, f_o
            ivectUI = [vars.ui.ei;vars.ui.eo];
            ivect = vars.e([bidTo,bidFrom]);
            ivectIDs = refs.ebvid([bidTo,bidFrom]);
            ovect = vars.f([bidTo,bidFrom]);
         elseif all(caus == [2,2]) % Inputs: f_i, f_o | Outputs: e_i, e_o
            ivectUI = [vars.ui.fi;vars.ui.fo];
            ivect = vars.f([bidTo,bidFrom]);
            ivectIDs = refs.fbvid([bidTo,bidFrom]);
            ovect = vars.e([bidTo,bidFrom]);
         elseif all(caus == [1,2]) % Inputs: e_i, f_o | Outputs: f_i, e_o
            ivectUI = [vars.ui.ei;vars.ui.fo];
            ivect = [vars.e(bidTo);vars.f(bidFrom)];
            ivectIDs = [refs.ebvid(bidTo),refs.fbvid(bidFrom)];
            ovect = [vars.f(bidTo);vars.e(bidFrom)];
         elseif all(caus == [2,1]) % Inputs: f_i, e_o | Outputs: e_i, f_o
            ivectUI = [vars.ui.fi;vars.ui.eo];
            ivect = [vars.f(bidTo);vars.e(bidFrom)];
            ivectIDs = [refs.fbvid(bidTo),refs.ebvid(bidFrom)];
            ovect = [vars.e(bidTo);vars.f(bidFrom)];
         else
            error('EID %u, ET %s, unexpected EC!!!',eid,etype);
         end;
         
         % Formulating the constituent equations
         % [e_i/f_i] = func_i(f_i/e_i,f_o/e_o)
         %  and
         % [e_o/f_o] = func_o(f_i/e_i,f_o/e_o)
         if strcmp(exprtype,'CC')
            srtsys.aes(aeid:(aeid+1)) = num2cell(ovect - expr * ivect);
         elseif strcmp(exprtype,'CE')
            srtsys.aes(aeid:(aeid+1)) = num2cell(ovect - ...
               subs(expr,ivectUI,ivect));
         elseif strcmp(exprtype,'NE')
            % Recording numeric indices
            srtsys.nvinpbvids(nid:(nid+1)) = getinpids(ivectIDs);
            % Constructing equations
            srtsys.nhs(nid:(nid+1)) = expr;
            srtsys.aes(aeid:(aeid+1)) = num2cell(ovect - vars.n(nid:(nid+1)));
            % Incrementing indices
            nid = nid + 2;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 2;
      
      % Modulated nonconservative 2-port coupler element
      % Assuming:
      % - exprtype = CMC, NMC, CME or NME
      % - Exactly two bonds are connected to the element
      % - One bond points to the element
      % - One bond points from the element
      case 'MR2'
         bidFrom = ebcm.from{eid};
         bidTo = ebcm.to{eid};
         exprid = refs.exprs(eid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         cambid = refs.causes(eid);
         caus = bgs.causes{cambid};
         
         % Checking causality
         if all(caus == [1,1]) % Inputs: e_i, e_o | Outputs: f_i, f_o
            ivectUI = [vars.ui.ei;vars.ui.eo];
            ivect = vars.e([bidTo,bidFrom]);
            ivectIDs = refs.ebvid([bidTo,bidFrom]);
            ovect = vars.f([bidTo,bidFrom]);
         elseif all(caus == [2,2]) % Inputs: f_i, f_o | Outputs: e_i, e_o
            ivectUI = [vars.ui.fi;vars.ui.fo];
            ivect = vars.f([bidTo,bidFrom]);
            ivectIDs = refs.fbvid([bidTo,bidFrom]);
            ovect = vars.e([bidTo,bidFrom]);
         elseif all(caus == [1,2]) % Inputs: e_i, f_o | Outputs: f_i, e_o
            ivectUI = [vars.ui.ei;vars.ui.fo];
            ivect = [vars.e(bidTo);vars.f(bidFrom)];
            ivectIDs = [refs.ebvid(bidTo),refs.fbvid(bidFrom)];
            ovect = [vars.f(bidTo);vars.e(bidFrom)];
         elseif all(caus == [2,1]) % Inputs: f_i, e_o | Outputs: e_i, f_o
            ivectUI = [vars.ui.fi;vars.ui.eo];
            ivect = [vars.f(bidTo);vars.e(bidFrom)];
            ivectIDs = [refs.fbvid(bidTo),refs.ebvid(bidFrom)];
            ovect = [vars.e(bidTo);vars.f(bidFrom)];
         else
            error('EID %u, ET %s, unexpected EC!!!',eid,etype);
         end;
         
         % Formulating the constituent equations
         % [e_i/f_i] = func_m_i(f_i/e_i,f_o/e_o,[t],[m])
         %  and
         % [e_o/f_o] = func_m_o(f_i/e_i,f_o/e_o,[t],[m])
         if strcmp(exprtype,'CMC')
            srtsys.aes(aeid:(aeid+1)) = num2cell(ovect - expr * ivect);
         elseif strcmp(exprtype,'NMC')
            % Retrieving modulating variable list
            modid = refs.mvars(eid);
            mvars = bgs.mvars{modid};
            % Identifying modulating variables
            [srtsys.nvmitids(nid:(nid+3)),srtsys.nvminpxids(nid:(nid+3)),...
               srtsys.nvminpbvids(nid:(nid+3))] = getvvids(mvars,vars,refs,...
               'vars2check','txb','vars2order',[false,true,true]);
            % Constructing equations
            srtsys.nhs(nid:(nid+3)) = expr.';
            nmat = reshape(vars.n(nid:(nid+3)),2,2).';
            srtsys.aes(aeid:(aeid+1)) = num2cell(ovect - nmat * ivect);
            % Incrementing indices
            nid = nid + 4;
         elseif strcmp(exprtype,'CME')
            srtsys.aes(aeid:(aeid+1)) = num2cell(ovect - ...
               subs(expr,ivectUI,ivect));
         elseif strcmp(exprtype,'NME')
            % Retrieving modulating variable list
            modid = refs.mvars(eid);
            mvars = bgs.mvars{modid};
            % Identifying modulating variables
            [srtsys.nvmitids(nid:(nid+1)),srtsys.nvminpxids(nid:(nid+1)),...
               srtsys.nvminpbvids(nid:(nid+1))] = getvvids(mvars,vars,refs,...
               'vars2check','txb','vars2order',[false,true,true]);
            % Recording numeric indices
            srtsys.nvinpbvids(nid:(nid+1)) = getinpids(ivectIDs);
            % Constructing equations
            srtsys.nhs(nid:(nid+1)) = expr;
            srtsys.aes(aeid:(aeid+1)) = num2cell(ovect - vars.n(nid:(nid+1)));
            % Incrementing indices
            nid = nid + 2;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + 2;
      
      % Nonconservative N-port coupler element
      % Assuming:
      % - exprtype = CC, CE or NE
      % - Exactly N bonds are connected to the element
      % - Expressions are given in the order of ports, not in the order of
      %   bonds
      case 'RN'
         exprid = refs.exprs(eid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         cambid = refs.causes(eid);
         caus = bgs.causes{cambid};
         pambid = refs.pnums(eid);
         pnum = bgs.epnums(pambid);
         bambid = refs.bcons(eid);
         bcon = bgs.ebcons{bambid};
         
         % Checking causality
         pidsEi = find(caus == 1);
         bidsEi = bcon(pidsEi);
         pidsFi = find(caus == 2);
         bidsFi = bcon(pidsFi);
         varsEFi = [vars.e(bidsEi);vars.f(bidsFi)];
         [~,rPIDs] = sort([pidsEi,pidsFi]);
         ivect = varsEFi(rPIDs);
         varsEFp = [vars.ep(pidsEi);vars.fp(pidsFi)];
         ivectUI = varsEFp(rPIDs);
         efBVIDs = [refs.ebvid(bidsEi),refs.fbvid(bidsFi)];
         ivectIDs = efBVIDs(rPIDs);
         varsEFo = [vars.f(bidsEi);vars.e(bidsFi)];
         ovect = varsEFo(rPIDs);
         
         % Formulating the constituent equations
         % [e_p1/f_p1] = func_p1(f_p1/e_p1,f_p2/e_p2,...)
         % [e_p2/f_p2] = func_p2(f_p1/e_p1,f_p2/e_p2,...)
         %  ...
         % [e_pN/f_pN] = func_pN(f_p1/e_p1,f_p2/e_p2,...)
         if strcmp(exprtype,'CC')
            srtsys.aes(aeid:(aeid+pnum-1)) = num2cell(ovect - expr * ivect);
         elseif strcmp(exprtype,'CE')
            srtsys.aes(aeid:(aeid+pnum-1)) = num2cell(ovect - ...
               subs(expr,ivectUI,ivect));
         elseif strcmp(exprtype,'NE')
            % Recording numeric indices
            srtsys.nvinpbvids(nid:(nid+pnum-1)) = getinpids(ivectIDs);
            % Constructing equations
            srtsys.nhs(nid:(nid+pnum-1)) = expr;
            srtsys.aes(aeid:(aeid+pnum-1)) = num2cell(ovect - ...
               vars.n(nid:(nid+pnum-1)));
            % Incrementing indices
            nid = nid + pnum;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + pnum;
      
      % Modulated nonconservative N-port coupler element
      % Assuming:
      % - exprtype = CMC, NMC, CME or NME
      % - Exactly N bonds are connected to the element
      % - Expressions are given in the order of ports, not in the order of
      %   bonds
      case 'MRN'
         exprid = refs.exprs(eid);
         exprtype = bgs.exprtypes{exprid};
         expr = bgs.exprs{exprid};
         cambid = refs.causes(eid);
         caus = bgs.causes{cambid};
         pambid = refs.pnums(eid);
         pnum = bgs.epnums(pambid);
         bambid = refs.bcons(eid);
         bcon = bgs.ebcons{bambid};
         
         % Checking causality
         pidsEi = find(caus == 1);
         bidsEi = bcon(pidsEi);
         pidsFi = find(caus == 2);
         bidsFi = bcon(pidsFi);
         varsEFi = [vars.e(bidsEi);vars.f(bidsFi)];
         [~,rPIDs] = sort([pidsEi,pidsFi]);
         ivect = varsEFi(rPIDs);
         varsEFp = [vars.ep(pidsEi);vars.fp(pidsFi)];
         ivectUI = varsEFp(rPIDs);
         efBVIDs = [refs.ebvid(bidsEi),refs.fbvid(bidsFi)];
         ivectIDs = efBVIDs(rPIDs);
         varsEFo = [vars.f(bidsEi);vars.e(bidsFi)];
         ovect = varsEFo(rPIDs);
         
         % Formulating the constituent equations
         % [e_p1/f_p1] = func_m_p1(f_p1/e_p1,f_p2/e_p2,...)
         % [e_p2/f_p2] = func_m_p2(f_p1/e_p1,f_p2/e_p2,...)
         %  ...
         % [e_pN/f_pN] = func_m_pN(f_p1/e_p1,f_p2/e_p2,...)
         if strcmp(exprtype,'CMC')
            srtsys.aes(aeid:(aeid+pnum-1)) = num2cell(ovect - expr * ivect);
         elseif strcmp(exprtype,'NMC')
            % Retrieving modulating variable list
            modid = refs.mvars(eid);
            mvars = bgs.mvars{modid};
            % Identifying modulating variables
            [srtsys.nvmitids(nid:(nid+pnum^2-1)),...
               srtsys.nvminpxids(nid:(nid+pnum^2-1)),...
               srtsys.nvminpbvids(nid:(nid+pnum^2-1))] = getvvids(mvars,vars,...
               refs,'vars2check','txb','vars2order',[false,true,true]);
            % Constructing equations
            srtsys.nhs(nid:(nid+pnum^2-1)) = expr.';
            nmat = reshape(vars.n(nid:(nid+pnum^2-1)),pnum,pnum).';
            srtsys.aes(aeid:(aeid+1)) = num2cell(ovect - nmat * ivect);
            % Incrementing indices
            nid = nid + pnum^2;
         elseif strcmp(exprtype,'CME')
            srtsys.aes(aeid:(aeid+pnum-1)) = num2cell(ovect - ...
               subs(expr,ivectUI,ivect));
         elseif strcmp(exprtype,'NME')
            % Retrieving modulating variable list
            modid = refs.mvars(eid);
            mvars = bgs.mvars{modid};
            % Identifying modulating variables
            [srtsys.nvmitids(nid:(nid+pnum-1)),...
               srtsys.nvminpxids(nid:(nid+pnum-1)),...
               srtsys.nvminpbvids(nid:(nid+pnum-1))] = getvvids(mvars,vars,...
               refs,'vars2check','txb','vars2order',[false,true,true]);
            % Recording numeric indices
            srtsys.nvinpbvids(nid:(nid+pnum-1)) = getinpids(ivectIDs);
            % Constructing equations
            srtsys.nhs(nid:(nid+pnum-1)) = expr;
            srtsys.aes(aeid:(aeid+pnum-1)) = num2cell(ovect - ...
               vars.n(nid:(nid+pnum-1)));
            % Incrementing indices
            nid = nid + pnum;
         else
            error('EID %u, ET %s, unexpected EXPRT!!!',eid,etype);
         end;
         aeid = aeid + pnum;
      
      % Error-catching case
      otherwise
         error('EID %u, unexpected ET!!!',eid);
   end;
end;

% Constructing the symbolic algebraic equations vector
srtsys.aes = sym(srtsys.aes);

% Converting numeric variable input tracking arrays to sparse matrices
srtsys.nvinpxids = inpids2mat(srtsys.nvinpxids,srtsys.xnum);
srtsys.nvinpbvids = inpids2mat(srtsys.nvinpbvids,srtsys.bvnum);
srtsys.nvminpxids = inpids2mat(srtsys.nvminpxids,srtsys.xnum);
srtsys.nvminpbvids = inpids2mat(srtsys.nvminpbvids,srtsys.bvnum);
srtsys.nvdepbvids = srtsys.nvinpbvids | srtsys.nvminpbvids;
srtsys.nvdepxids = srtsys.nvinpxids | srtsys.nvminpxids;
end

%% Input parser construction function

function ip = buildIP()
% The purpose of this function is to construct an input parser object for the
% formulateAEqs function.
%
% INPUTS/OUTPUTS
% ip - Input parser object for the formulateAEqs function.

ip = inputParser;
ip.StructExpand = false;
ip.addRequired('bgs',@isstruct);
ip.addRequired('refs',@isstruct);
ip.addRequired('vars',@isstruct);
ip.addRequired('ebcm',@isstruct);
end

%% NV input tracking array conversion function

function inpidsMat = inpids2mat(inpids,trkvarnum)
% The purpose of this function is to convert an NV input tracking array from an
% array of 1-/2-field data structures to a sparse boolean/double input tracking
% array.
%
% INPUTS/OUTPUTS
% inpids - NV input tracking array, a vector of 1-/2-field input tracking data
%          structures. inpids can be tracking XID inputs (regular and modulated)
%          and BVID inputs (regular and modulated). It cannot be tracking time
%          modulation, as that tracking array is already a boolean vector. This
%          function checks whether or not the array is ordered (all of them
%          should be though).
% trkvarnum - Number of trackable variables of inpids' tracked type in the BGS.
% inpidsBool - NV input tracking array, a sparse boolean/double array.

% Constructing tracked VIDs vector
trkvids = double([inpids.indices]);
vnum = length(trkvids);

% Constructing NVIDs vector
nvids = zeros(1,vnum);
vid = 1; % ID in variable indices vector
for i = 1:length(inpids)
   L = length(inpids(i).indices);
   if L > 0
      nvids(vid:(vid+L-1)) = i;
      vid = vid + L;
   end;
end;

% Constructing order indices vector
if isfield(inpids,'order')
   oids = double([inpids.order]);
end;

% Converting NV input tracking array
if isfield(inpids,'order')
   inpidsMat = sparse(nvids,trkvids,oids,length(inpids),double(trkvarnum));
else
   inpidsMat = sparse(nvids,trkvids,true,length(inpids),double(trkvarnum));
end;
end

%% NV input 2-field data structures construction function

function inpids = getinpids(indices,iforder,order)
% The purpose of this function is to construct a 2-field data structure (or
% 1-field, if the inputs are not ordered) with .indices (and optionally, .order)
% field(s). By default, the function considers the input indices vector to be in
% the appropriate (monotonically increasing, [1,2,...]) order. It sorts the
% indices if necessary, and appropriately rearranges the order vector.
%
% INPUTS/OUTPUTS
% indices - Horizontal full uint16 vector of input variable IDs. REQUIRED.
% iforder - Whether or not to output the .order field. OPTIONAL. By default, the
%           field is included.
% order - Horizontal full uint8 in-input order vector. OPTIONAL. By default, the
%         field is included, and assumes that indices are in the appropriate
%         order as input.
% inpids - 1-/2-field input tracking data structure. It has the following
%          fields:
%  .indices - Horizontal full uint16 vector of variable IDs, sorted in ascending
%             order.
%  .order - Horizontal full uint8 vector of variable positions in the input
%           vector, only present if the input is ordered.

if nargin == 1
   iforder = true;
   order = uint8(1):uint8(length(indices));
elseif nargin == 2
   if iforder
      order = uint8(1):uint8(length(indices));
   end;
elseif nargin == 3
   order = uint8(order);
else
   error('Invalid number of inputs!!!');
end;
indices = uint16(indices);

[inpids.indices,ind] = sort(indices);
if iforder
   inpids.order = order(ind);
end;
end