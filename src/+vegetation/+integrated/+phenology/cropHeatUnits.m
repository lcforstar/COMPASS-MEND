function [CHU,HU,HUI,DRDS,RDS,HUF,dHUF] = cropHeatUnits( ...
        Tmax,Tmin,dates,sowingDoy,baseT,PHU,F1,L1,F2,L2)
%CROPHEATUNITS Recurring calendar-based EPIC heat-unit forcing auxiliary.
% Temperatures [degC, nday x 1], dates [datetime, nday x 1], sowingDoy [DOY],
% baseT [degC], PHU [degC d] and existing EPIC dimensionless S-curve points.
% Output arrays are nday x 1; CHU/HU [degC d], others dimensionless. Each
% configured sowing DOY starts a crop season; January 1 never resets it.
% Maturity saturates this empirical curve until the next sowing event. The
% curve diagnoses LAI potential, never manages physical harvest or C/N pools.
% Equations are the existing vegetation.phenology.CalculateHeatUnit EPIC
% closure. Calendar adaptation corrects the DOY/row-index ambiguity; dHUF is
% today's nonnegative increment, with zero at each crop-season reset.
% Editor: Wenjie Liu; Date: 2026-09-30.
    n=numel(Tmax);Tmax=Tmax(:);Tmin=Tmin(:);dates=dates(:);
    assert(numel(Tmin)==n && numel(dates)==n && isdatetime(dates), ...
        'cropHeatUnits:shape','Daily dates and extrema must have equal lengths.');
    assert(isfinite(PHU)&&PHU>0 && sowingDoy>=1&&sowingDoy<=366, ...
        'cropHeatUnits:parameters','Positive PHU and valid sowing DOY are required.');
    l2=(log(F1*(1-L1)/L1)-log(F2*(1-L2)/L2))/(F2-F1);
    l1=log(F1*(1-L1)/L1)+l2*F1;
    CHU=zeros(n,1);HU=CHU;HUI=CHU;DRDS=CHU;RDS=CHU;HUF=CHU;dHUF=CHU;
    active=false;sumHu=0;previousHuf=0;mature=false;
    for d=1:n
        planting=day(dates(d),'dayofyear')==sowingDoy;
        if planting,active=true;sumHu=0;previousHuf=0;mature=false;end
        if ~active,continue;end
        if ~mature
            HU(d)=max((Tmax(d)+Tmin(d))/2-baseT,0);
            sumHu=sumHu+HU(d);
            DRDS(d)=HU(d)/PHU;
        end
        CHU(d)=sumHu;HUI(d)=sumHu/PHU;RDS(d)=HUI(d);
        HUF(d)=HUI(d)/(HUI(d)+exp(l1-l2*HUI(d)));
        if ~planting,dHUF(d)=max(0,HUF(d)-previousHuf);end
        previousHuf=HUF(d);
        mature=mature||HUI(d)>=1;
    end
end
