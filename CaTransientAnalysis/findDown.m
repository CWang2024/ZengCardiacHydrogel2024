function [OutTime, OutTrace, OutDouble] = findDown(Trace, Times, ratio, DownProm, fps)

%IMPORTNAT: Trace and Time here can only take rising phase input;

PreciseAmp=Trace(1) - ratio*DownProm;
DownInd1 = find(Trace >= PreciseAmp,1,'last');

%DownInd2 = find(Trace <= PreciseAmp,1,'first');
DownInd2 = find(Trace <= PreciseAmp);
DownInd2 = find(DownInd1 <=DownInd2, 1, 'first');
k=Trace(DownInd2)-Trace(DownInd1);
b=Trace(DownInd1)-k*DownInd1;
OutDouble = (PreciseAmp-b)/k;
OutTime = Times(DownInd1)+(OutDouble-DownInd1)/fps;
OutTrace = PreciseAmp;
end

