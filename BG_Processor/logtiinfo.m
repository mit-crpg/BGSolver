%% logtiinfo
% Time integrator information logging function.
%
% logtiinfo(fids,runopts,tioptions,ti_info)
%
% This function logs to specification (using log2spec) the known (if any)
% information about the specified time integrator and the integration options.
% It relies on the time integrator database, stored in a .mat file in the
% mat_Files subpackage.
%
% Package:    BGSolver v1.03
% Subpackage: BG_Processor
% Date:       November 8, 2012
% Author:     Eugeny Sosnovsky
%             esos@mit.edu

%% Time integrator information logging function

function logtiinfo(fids,runopts,tioptions,ti_info)
% The purpose of this function is to log2spec the known (if any) information
% about the time integrator and the integration options.
%
% INPUTS/OUTPUTS
% fids - File IDs data structure.
% runopts - Runtime and Recording Options data structure.
% tioptions - Time Integrator options data structure.
% ti_info - Time Integrator database data structure. Loaded from ti_info.mat.

log2spec(sprintf('Time Integrator: %s',runopts.tit),1,fids,runopts);

if isfield(ti_info,runopts.tit)
   ti = ti_info.(runopts.tit);
   log2spec(sprintf('%s found in the BGSolver Time Integrator database.',...
      runopts.tit),1,fids,runopts);
   log2spec(sprintf('Source:             %s',ti.src),2,fids,runopts);
   log2spec(sprintf('Stiffness:          %s',ti.stiff),2,fids,runopts);
   log2spec(sprintf('Order of accuracy:  %s',ti.order),2,fids,runopts);
   log2spec(sprintf('Explicitness:       %s',ti.explicit),2,fids,runopts);
   log2spec(sprintf('Method type:        %s',ti.method),2,fids,runopts);
   log2spec(sprintf('Scaling:            %s',ti.scaling),2,fids,runopts);
   log2spec(sprintf('Vectorization:      %s',ti.vector),2,fids,runopts);
   log2spec(sprintf('Operator splitting: %s',ti.opsplit),2,fids,runopts);
   log2spec(sprintf('Time step:          %s',ti.dt),2,fids,runopts);
   log2spec(sprintf('Runtime recording:  %s',ti.rr),2,fids,runopts);
   log2spec(sprintf('Jumpstarter types:  %s',ti.jstype),2,fids,runopts);
   log2spec(sprintf('Residual scaling:   %s',ti.rscal),2,fids,runopts);
   log2spec(sprintf('X vector scaling:   %s',ti.xscal),2,fids,runopts);
   log2spec(sprintf('X vector centering: %s',ti.xcent),2,fids,runopts);
   log2spec(sprintf('Preconditioner:     %s',ti.pc),2,fids,runopts);
   log2spec(sprintf('Predictor:          %s',ti.igpdct),2,fids,runopts);
   log2spec(sprintf('Jacobian:           %s',ti.jacob),2,fids,runopts);
   log2spec(sprintf('Additional info:    %s',ti.addinfo),2,fids,runopts);
   if tioptions.phxsys.phxsplit > 0
      if ti.ifopsplit
         log2spec(['Split physics option enabled, time integrator does ',...
            'support physics splitting.'],1,fids,runopts);
         if tioptions.phxsys.phxsplit == 1
            cpl = 'Triangular';
         elseif tioptions.phxsys.phxsplit == 2
            cpl = 'Diagonal';
         end;
         log2spec(sprintf('%s coupling will be used.',cpl),1,fids,runopts);
      else
         log2spec(['Split physics option enabled, but time integrator does ',...
            'not support physics splitting.'],1,fids,runopts);
         log2spec('Full coupling will be forced instead.',1,fids,runopts);
      end;
   else
      log2spec('Fully coupled option enabled.',1,fids,runopts);
   end;
   if ti.ifvectreq
      log2spec('Time integrator does require vectorization.',2,fids,runopts);
   else
      log2spec('Time integrator does not require vectorization.',2,fids,...
         runopts);
   end;
   if ti.Njs == 0
      log2spec('Time integrator does not require a jumpstart.',2,fids,runopts);
   else
      log2spec('Time integrator does require a jumpstart.',2,fids,runopts);
      if ti.pts
         log2spec(['Time integrator can use pre-initial time steps if ',...
            'provided.']);
      else
         log2spec(['Time integrator cannot use pre-initial time steps if ',...
            'provided.'],2,fids,runopts);
         log2spec('Zero-derivative jumpstart will be used if necessary.',...
            2,fids,runopts);
      end;
      if runopts.jstype == -1
         log2spec(sprintf('Default %s jumpstart enabled.',ti.defjs),...
            2,fids,runopts);
      elseif runopts.jstype == 0
         log2spec('Zero-derivative jumpstart enabled.',2,fids,runopts);
      elseif runopts.jstype == 1
         log2spec('Ministep-based jumpstart enabled.',2,fids,runopts);
      end;
   end;
else
   log2spec(sprintf('%s is not in the BGSolver Time Integrator database.',...
      runopts.tit),1,fids,runopts);
   if tioptions.phxsys.phxsplit > 0
      log2spec(['Split physics option enabled, but cannot confirm split ',...
         'physics support.'],1,fids,runopts);
      log2spec('The integrator may force full coupling.',1,fids,runopts);
   else
      log2spec('Fully coupled option enabled.',1,fids,runopts);
   end;
   log2spec('Vectorization enabled as a precaution.',2,fids,runopts);
   log2spec('Time integrator may require a jumpstart.',2,fids,runopts);
   log2spec('Zero-derivative jumpstart will be used if necessary.',...
      2,fids,runopts);
end;

if runopts.tst == 0
   log2spec('Adaptive time stepping enabled.',1,fids,runopts);
elseif runopts.tst == 1
   log2spec('Fixed time stepping enabled.',1,fids,runopts);
   log2spec(sprintf('Fixed step size %g.',tioptions.tdisc.dt),1,fids,runopts);
elseif runopts.tst == 2
   log2spec(['Fixed time stepping with sequential adaptive integrations ',...
      'enabled.'],1,fids,runopts);
   log2spec(sprintf('Fixed step size %g.',tioptions.tdisc.dt),1,fids,runopts);
end;
end