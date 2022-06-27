finewake = 1;
% define wakefields
if finewake
  decm=1;
  bw=2.5e-3;
else
  decm=10;
  bw=0.1;
end
% -- S-band
load /Users/cemma/Documents/Work/FACET-II/Lucretia_sims/facet2-lattice/Lucretia/data/srwf_long_sband.mat wf wfz
WF.ZSR(1).z=wfz(1:decm:end); WF.ZSR(1).K=wf(1:decm:end); WF.ZSR(1).BinWidth=bw;
load /Users/cemma/Documents/Work/FACET-II/Lucretia_sims/facet2-lattice/Lucretia/data/srwf_trans_sband.mat wf wfz
WF.TSR(1).z=wfz(1:decm:end); WF.TSR(1).K=wf(1:decm:end); WF.TSR(1).BinWidth=bw;
% -- X-band
load /Users/cemma/Documents/Work/FACET-II/Lucretia_sims/facet2-lattice/Lucretia/data/srwf_long_xband.mat wf wfz
WF.ZSR(2).z=wfz; WF.ZSR(2).K=wf; WF.ZSR(2).BinWidth=bw;
load /Users/cemma/Documents/Work/FACET-II/Lucretia_sims/facet2-lattice/Lucretia/data/srwf_trans_xband.mat wf wfz
WF.TSR(2).z=wfz; WF.TSR(2).K=wf; WF.TSR(2).BinWidth=bw;

% assign wakefields
id=findcells(BEAMLINE,'Class','LCAV',1,findcells(BEAMLINE,'Name','BEGBC20'));
for n=1:length(id)
  if isempty(regexp(BEAMLINE{id(n)}.Name,'^L1XF','once')) % s-band
    BEAMLINE{id(n)}.Wakes=[1,1];
  else % x-band (comment out for now)
%     if phaseInp(4)==0
%       BEAMLINE{id(n)}.Wakes=[0,0];
%     else
%       BEAMLINE{id(n)}.Wakes=[2,2];
%     end
%     BEAMLINE{id(n)}.Phase=180;
%     BEAMLINE{id(n)}.Volt=phaseInp(4);
%     BEAMLINE{id(n)}.Egain=phaseInp(4)*cos(BEAMLINE{id(n)}.Phase);
  end
end

% activate longitudinal wakes
SetTrackFlags('SRWF_Z',1,1,findcells(BEAMLINE,'Name','BEGBC20'));

% activate transverse wakesTwiss
SetTrackFlags('SRWF_T',0,1,findcells(BEAMLINE,'Name','BEGBC20'));