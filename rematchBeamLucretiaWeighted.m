function [beamStructOutput]=rematchBeamLucretiaWeighted(mbetax,mbetay,malphax,malphay,beamStruct)
beamStructOutput = beamStruct;
ps = beamStruct.Bunch.x;
w = (beamStruct.Bunch.Q/sum(beamStruct.Bunch.Q))*length(beamStruct.Bunch.Q);% Weighted charge
%[x xp y yp z E(GeV)] are coordinates of ebeamin
    x=ps(1,:)-mean(w.*ps(1,:));
    xp=ps(2,:)-mean(w.*ps(2,:));
    y=ps(3,:)-mean(w.*ps(3,:));
    yp=ps(4,:)-mean(w.*ps(4,:));    
    gamma = ps(6,:)/0.511e-3;
    gavg = mean(w.*gamma);
    % Input beam Twiss parameters
    emitx=sqrt(mean(w.*x.^2).*mean(w.*xp.^2)-mean(w.*x.*xp).^2).*gavg       
    emity=sqrt(mean(w.*y.^2).*mean(w.*yp.^2)-mean(w.*y.*yp).^2).*gavg
    betax=mean(w.*x.*x).*gavg./emitx
    betay=mean(w.*y.*y).*gavg./emity
    alphax=-mean(w.*x.*xp).*gavg./emitx
    alphay=-mean(w.*y.*yp).*gavg./emity

    %remove old correlation
    xp=xp+alphax.*x./betax;
    yp=yp+alphay.*y./betay;
    % From here on down the x,xp,y,yp co-ordinates are the new re-matched coordinates
    %scale the beam
    x=x.*sqrt(mbetax./betax);
    y=y.*sqrt(mbetay./betay);
    xp=xp.*sqrt(betax./mbetax);
    yp=yp.*sqrt(betay./mbetay);

    %add new correlation
    xp=xp-malphax.*x./mbetax;
    yp=yp-malphay.*y./mbetay;

    % Output beam Twiss parameters
    emitxnew=sqrt(mean(w.*x.^2).*mean(w.*xp.^2)-mean(w.*x.*xp).^2).*gavg
    emitynew=sqrt(mean(w.*y.^2).*mean(w.*yp.^2)-mean(w.*y.*yp).^2).*gavg
    betaxnew=mean(w.*x.*x).*gavg./emitxnew
    betaynew=mean(w.*y.*y).*gavg./emitynew
    alphaxnew=-mean(w.*x.*xp).*gavg./emitxnew
    alphaynew=-mean(w.*y.*yp).*gavg./emitynew
    
    % Output beam structure
    beamStructOutput.Bunch.x(1,:) = x;
    beamStructOutput.Bunch.x(2,:) = xp;
    beamStructOutput.Bunch.x(3,:) = y;
    beamStructOutput.Bunch.x(4,:) = yp;
    