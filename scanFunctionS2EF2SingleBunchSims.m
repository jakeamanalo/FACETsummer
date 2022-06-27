function [data,BEAMLINE] = scanFunctionS2EF2SingleBunchSims(scanparams)
%SCANFUNCTIONNOTCHCOLLSIMS Summary of this function goes here
%   Detailed explanation goes here
addpath(genpath('/Users/cemma/Documents/Work/FACET-II/Laser_work/Laser_heater/LaserHeaterSimulations'))
addpath(genpath('/Users/cemma/Documents/Work/FACET-II/Lucretia_sims/facet2-lattice'))
% load the lattice - this loads the entire 1 km lattice
load '/Users/cemma/Documents/Work/FACET-II/Lucretia_sims/facet2-lattice/Lucretia/models/FACET2e/FACET2e'
% Load the electron beam for the 2nC single bunch case
%load '/Users/cemma/Documents/Work/FACET-II/Lucretia_sims/facet2-lattice/Lucretia/beams/FACET2e_1M.mat'
Beam = importdata('matchedBeamAtENDINJ.mat');
ISRon=0; % 1= Incoherent Sync. Rad. on, 0= ISR off
CSRon=0; % 1= Coherent Sync. Rad. on, 0= CSR off
LSCon=0; % 1= Longitudinal space charge on, 0= LSC off
% Decimate beam?
decBeam=100; % use integer>1 
             
% -- Tracking options
% Define indices
istart=findcells(BEAMLINE,'Name','ENDINJ');
iLH=findcells(BEAMLINE,'Name','HTRUNDF');
iL1A=findcells(BEAMLINE,'Name','BEGL1F');
BC11END=findcells(BEAMLINE,'Name','BC11CEND'); % end of BC11
begBC14 =findcells(BEAMLINE,'Name','BEGBC14_1');
endBC14 =findcells(BEAMLINE,'Name','ENDBC14_2');
iBC1B=findcells(BEAMLINE,'Name','ENDBC11_2'); % start of L2
iBC2B=findcells(BEAMLINE,'Name','BEGL3F_1'); % start of L3
iL3B=findcells(BEAMLINE,'Name','ENDL3F_2'); % end of L3
BEGBC20 =findcells(BEAMLINE,'Name','BEGBC20');
ENDBC20 =findcells(BEAMLINE,'Name','ENDBC20');
iFF=findcells(BEAMLINE,'Name','ENDBC20'); % start of FFS
QM10631=findcells(BEAMLINE,'Name','QM10631'); % start of FFS
QM10631 = QM10631(1);
% -- Add wakefields
BEAMLINE{194}.Wakes = [0 0];% This structure had wakes set to [4 0] which made the tracking crash
loadFacetWakefields
% Collective effects
%-- Sync. Radiation
eleSR=findcells(BEAMLINE,'Class','SBEN');
for iele=eleSR
  BEAMLINE{iele}.TrackFlag.SynRad=2*ISRon;
  BEAMLINE{iele}.TrackFlag.CSR=-1*CSRon;
  BEAMLINE{iele}.TrackFlag.CSR_SmoothFactor=0;
  BEAMLINE{iele}.TrackFlag.CSR_DriftSplit=25;
  BEAMLINE{iele}.TrackFlag.Split=25;
  % Turn off 2d CSR
        BEAMLINE{iele}.TrackFlag.CSR_2D = 0;
        BEAMLINE{iele}.TrackFlag.CSR_USEGPU = 0;
end

% -- Space charge
for iele=findcells(BEAMLINE,'TrackFlag')
  BEAMLINE{iele}.TrackFlag.LSC=LSCon;
  BEAMLINE{iele}.TrackFlag.LSC_storeData=0;
  % Set NBPM on LCAV elements to ensure 0.1m drift sections for
  % application of LSC
  if strcmp(BEAMLINE{iele}.Class,'LCAV')
    BEAMLINE{iele}.NBPM=LSCon*BEAMLINE{iele}.L/0.1;
    BEAMLINE{iele}.GetSBPMData=LSCon;
    BEAMLINE{iele}.GetInstData=LSCon;
  end
end

% Set initial beam parameters
Initial.Q=Initial.Q*(1+scanparams.qi); % Initial charge offset
Beam.Bunch.x(1,:)=Beam.Bunch.x(1,:)+scanparams.dx;% Initial offset in x
Beam.Bunch.x(3,:)=Beam.Bunch.x(3,:)+scanparams.dy;% Initial offset in y

% --- Set Linac Phases / degrees off-crest and amplitudes
[BEAMLINE]=changeLinacPhases(BEAMLINE,scanparams);

% Decimate beam
if decBeam>1
  decBeam=floor(decBeam);
  Beam.Bunch.x=Beam.Bunch.x(:,1:decBeam:end);
  Beam.Bunch.stop=Beam.Bunch.stop(1:decBeam:end);
  Beam.Bunch.Q=Beam.Bunch.Q(1:decBeam:end).*decBeam;
  disp('Beam decimated')
end   
%% Track through the linac until BC20BEG
% Track to end of each bending section and re-center beam in each case to
% take care of phase-slip with respect to RF and orbit excursions due to SR
% energy losses. Also remove linear dispersion induced by CSR effects.
% In reality this is done with RF phasing, orbit feedbacks and beam tuning.
% Track through the laser heater
[~,beam]=TrackThru(istart,iLH,Beam,1,1,0); % Track to Laser heater and add energy spread
LH = 0;
beam.Bunch.x(6,~beam.Bunch.stop)=beam.Bunch.x(6,~beam.Bunch.stop)+randn(size(beam.Bunch.x(6,~beam.Bunch.stop))).*LH.*1e-6;
% figure;
% subplot(2,1,1)
% scatter(Beam.Bunch.x(5,~beam.Bunch.stop),Beam.Bunch.x(6,~beam.Bunch.stop),'.');
% subplot(2,1,2)
% scatter(beam.Bunch.x(5,~beam.Bunch.stop),beam.Bunch.x(6,~beam.Bunch.stop),'.r');
%Continue tracking form the Laser Heater to the end of BC20
ele=[iLH iL1A iBC1B iBC2B BEGBC20 ENDBC20];
for iele=1:length(ele)-1
  [~,beam]=TrackThru(ele(iele),ele(iele+1),beam,1,1);
  disp(BEAMLINE{ele(iele)}.Name)
  for ii=1:5
    if ii<5
      coef=polyfit(beam.Bunch.x(6,~beam.Bunch.stop),beam.Bunch.x(ii,~beam.Bunch.stop),1);
      coef(2)=0;
      beam.Bunch.x(ii,~beam.Bunch.stop)=beam.Bunch.x(ii,~beam.Bunch.stop)-polyval(coef,beam.Bunch.x(6,~beam.Bunch.stop));
    end
    beam.Bunch.x(ii,~beam.Bunch.stop)=beam.Bunch.x(ii,~beam.Bunch.stop)-mean(beam.Bunch.x(ii,~beam.Bunch.stop));
  end
  beamImage(beam);
  [nx,ny,~] = GetNEmitFromBeam( beam, 1,1 );
  disp(sprintf(['Emit at ',BEAMLINE{ele(iele+1)}.Name,' = ',num2str(nx),', ', num2str(ny)]))
end
%% Plot the output LPS at the exit of BC20
        ebins = [9.5:0.01:10.2];% GeV % The resolution of ~0.5 MeV/pix comes from the pixel size/dispersion at TCAV
        zbins = ([-150:0.5:100]*1e-6);% m % The resolution of ~0.25 um/pix comes from the pixel size/streak in um/um
        [~,~,~,data.I,data.Eprof]=MakeBeamLPS(beam,zbins,ebins,1);        
end

