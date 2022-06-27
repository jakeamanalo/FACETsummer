% Simple Lucretia Tutorial
% Simulates a single RF accelerating cavity followed by a 4 dipole bunch
% bunch compressor
%%%%%%%%%%%%%%%%%%%% PARAMETERS YOU CAN PLAY WITH %%%%%%%%%%%%%%%%%%%%%%%%%
initialEnergy = .5;                 % E beam energy [GeV]
initialBunchLength = 10e-6;         % E beam bunch length [m]
bendAngle = 2.5e-2;                 % Bend angle for each dipole in chicane [rad]
lbend = 2;                         % Dipole length for each dipole in chicane [m]
ldrift = 2;                         % Drift length between dipoles in chicane [m]
rfphaseAngle = -55;                 % RF Cavity Phase [V]
rfVoltage = 3e3;                    % RF Cavity Voltage [V]
ISRon = 0;CSRon = 0;LSCon = 0;      % Special Physics effects 0= off 1=on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global BEAMLINE
%% Creates the Initial Lucretia beam and Initial Condition structure
[beam,Initial]=generateInitialLucretiaBeam(initialEnergy,initialBunchLength);
%% Creates the beamline: 1 RF cavity 3 quads, 4 bends, 8 drifts
[BEAMLINE] = generateLucretiaBeamline(Initial,bendAngle,rfVoltage,rfphaseAngle,lbend,ldrift);
%% Add 'special' physics effects to the simulation
% -- CSR
    eleSR=findcells(BEAMLINE,'Class','SBEN');
    for iele=eleSR
      BEAMLINE{iele}.TrackFlag.SynRad=2*ISRon;
      BEAMLINE{iele}.TrackFlag.CSR=-1*CSRon;
      BEAMLINE{iele}.TrackFlag.CSR_SmoothFactor=1;
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
  BEAMLINE{iele}.TrackFlag.ZMotion =1;
    end  
%% Run the simulation stopping after the accelerating cavity and at the end     
tic
    T = Track(beam);
    T.startInd=1; T.finishInd=2; T.trackThru();
    beamAfterAcceleratingCavity=T.beamOut;   
    beamImage(beamAfterAcceleratingCavity);
    T = Track(beamAfterAcceleratingCavity);
    T.startInd = 2;T.finishInd=length(BEAMLINE); T.trackThru();
    beamAtExit=T.beamOut;   
    beamImage(beamAtExit);   
toc

function [beam,Initial]= generateInitialLucretiaBeam(initialEnergy,initialBunchLength)
Initial = InitCondStruc( );
% Set the energy and bunch length
    Initial.Momentum = initialEnergy;% GeV
    Initial.sigz = initialBunchLength;% m
% Set the beam Twiss parameters
    Initial.x.Twiss.beta = 10;  Initial.x.Twiss.alpha = 0.8;
    Initial.y.Twiss.beta = 10;  Initial.y.Twiss.alpha = 0.8;
    Initial.x.NEmit = 1e-6;     Initial.y.NEmit = 1e-6;
% Set intial energy spread and a chirp 
    Initial.SigPUncorrel = 1e-3*Initial.Momentum;
    Initial.PZCorrel = 0*Initial.Momentum;% 0 initial chirp
% Set the Charge
    Initial.Q = 2e-9;
% Make a 6D Gaussian beam with the given Initial conditions
    beam = MakeBeam6DGauss( Initial, 1e4, 5, 1 );
    beamImage(beam)
end

function [BEAMLINE] = generateLucretiaBeamline(Initial,bendAngle,rfVoltage,rfphaseAngle,lbend,ldrift)
% One RF Accelerating Cavity
    L = 2;    BEAMLINE{1} =RFStruc(L, rfVoltage, 0, 2856, 0, 0, 0, 5 , 'RF0' );   
    BEAMLINE{1}.Phase = rfphaseAngle;
    Efinal = Initial.Momentum + 1e-3*rfVoltage*cosd(rfphaseAngle);
    brho = Efinal/0.29979;% A constant
%Drift 1
    L = 0.5;    BEAMLINE{3} = DrifStruc(L, 'GAP');
%Quad 1
    L = 0.15;    K = 0/L/brho;   Aper = 1e-2;
    BEAMLINE{2} = QuadStruc(L, L*K*brho,0,Aper,'Q1');
%Drift 2
    L = 0.5;    BEAMLINE{3} = DrifStruc(L, 'GAP');
%Quad 2
    L = 0.15;    K = 0/L/brho;   Aper = 1e-2;
    BEAMLINE{4} = QuadStruc(L, L*K*brho,0,Aper,'Q2');
%Drift 3
    L = 0.5;    BEAMLINE{5} = DrifStruc(L, 'GAP');
%Quad 3
    L = 0.15;    K = 0/L/brho;   Aper = 1e-2;
    BEAMLINE{6} = QuadStruc(L, L*K*brho,0,Aper,'Q3');
%Drift 4
    L = 0.4;    BEAMLINE{7} = DrifStruc(L, 'DF2');    
%Bend 1
    L = lbend; Angle = bendAngle; B = Angle*brho;%[T.m]  
    hgap = 2.16e-2;
    BEAMLINE{8} = SBendStruc( L, B, Angle, 0, 0, hgap, 0, 0, 'B1' );
%Drift 5
    L = ldrift;    BEAMLINE{9} = DrifStruc(L, 'CSD1');  
%Bend 2
    L = lbend; Angle = -bendAngle; B = Angle*brho;%[T.m]  
    hgap = 2.16e-2;
    BEAMLINE{10} = SBendStruc( L, B, Angle, 0, 0, hgap, 0, 0, 'B2' );   
%Drift 6
    L = 0.25;    BEAMLINE{11} = DrifStruc(L, 'CSD3'); 
%Bend 3
    L = lbend; Angle = -bendAngle; B = Angle*brho;%[T.m]  
    hgap = 2.16e-2;
    BEAMLINE{12} = SBendStruc( L, B, Angle, 0, 0, hgap, 0, 0, 'B3' ); 
%Drift 7
    L = ldrift;    BEAMLINE{13} = DrifStruc(L, 'CSD1'); 
%Bend 4
    L = lbend; Angle = bendAngle; B = Angle*brho;%[T.m]  
    hgap = 2.16e-2;
    BEAMLINE{14} = SBendStruc( L, B, Angle, 0, 0, hgap, 0, 0, 'B4' );    
%Drift 8
    L = 0.25;    BEAMLINE{15} = DrifStruc(L, 'CSD3');     
% Set Central Energy of Beamline
    
    BEAMLINE{1}.P = Initial.Momentum;
    for n=1:length(BEAMLINE); BEAMLINE{n}.P = Efinal;end

% Set Positions of magnets
    SetSPositions(1,length(BEAMLINE),0)
% Look at the beamline and the Twiss Parameters
%    TwissPlot(1,length(BEAMLINE),Initial ,[-1 -1 0]);

end