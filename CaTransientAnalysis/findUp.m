function [OutTime, OutTrace, OutDouble] = findUp(Trace, Times, ratio, UpProm, fps)

%IMPORTNAT: Trace and Time here can only take rising phase input;

PreciseAmp=Trace(1) + ratio.*UpProm;
UpInd1 = find(Trace <= PreciseAmp);
UpInd2 = find(Trace >= PreciseAmp,1,'first');
UpInd1 = find(UpInd1 <= UpInd2, 1, 'last');
k=Trace(UpInd2)-Trace(UpInd1);
b=Trace(UpInd1)-k*UpInd1;
OutDouble = (PreciseAmp-b)/k;
OutTime = Times(UpInd1)+(OutDouble-UpInd1)/fps;
OutTrace = PreciseAmp;
end

