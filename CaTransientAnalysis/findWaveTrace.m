function waveTrace = findWaveTrace(MinInd1, MinInd2, IntData, Times)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
            disp("Entered");
            TraceInd = MinInd1:MinInd2;
            TraceTimes = Times(TraceInd);
            Trace = IntData(TraceInd);
            %OriginalTrace=OriginalIntData(TraceInd);
    
            MaxInt = max(Trace);
            MaxIntInd = find(Trace == MaxInt,1,'first');
            
            UpProm = MaxInt - Trace(1);
            DownProm = MaxInt - Trace(end);
            StartInd=1;
            DiffTrace=diff(Trace);
            for j=1:MaxIntInd-3
                %if (Trace(j+4)-Trace(j))>0.15*UpProm && Trace(j+1)-Trace(j)>0 && Trace(j+2)-Trace(j+1)>0
                if all(DiffTrace(j:j+2)>0)&& Trace(j+2)-Trace(j)>0.2*UpProm
                   StartInd=j;
                   break;
                end
            end
            

            EndInd=length(Trace);
            
            Down70Ind = find(Trace >= Trace(end) + 0.3*DownProm,1,'last')
            length(Trace)-2
            for n=Down70Ind+2:length(Trace)-3
                disp("haahahahahaha");
                if (Trace(n+2)-Trace(n)>0.3*(Trace(n)-Trace(n-2)))&&(DiffTrace(n)>-0.15)&&DiffTrace(n)<DiffTrace(n+1)
                    EndInd=n;
                    break;
                end
            end    
            %}

            %{
           
            if isempty(waveBaseInd)
                waveBaseInd= [waveBaseInd , (MinInd(t-1)+StartInd-1: MinInd(t-1)+EndInd-1)];
            else
                if waveBaseInd(end)== MinInd + StartInd -1
                    waveBaseInd= [waveBaseInd , (MinInd(t-1)+StartInd: MinInd(t-1)+EndInd-1)];
                else
                    waveBaseInd= [waveBaseInd , (MinInd(t-1)+StartInd-1: MinInd(t-1)+EndInd-1)];
                end
            end
            %}
                 
            %MaxInd=MaxIntInd-StartInd+1;
            %MinInt = MinValues(t-1);
            %MinIntInd = find(Trace == MinInt,1,'first');
            waveTrace=struct();
            waveTrace.Trace=Trace(StartInd: EndInd);
            waveTrace.TraceUp=Trace(StartInd:MaxIntInd);
            waveTrace.TraceDown=Trace(MaxIntInd:EndInd);
            waveTrace.Time=TraceTimes(StartInd: EndInd);
            waveTrace.TimeUp=TraceTimes(StartInd:MaxIntInd);
            waveTrace.TimeDown=TraceTimes(MaxIntInd:EndInd);
            waveTrace.Startloc=MinInd1+StartInd-1;
            waveTrace.Endloc=MinInd1+EndInd-1;
            waveTrace.Maxloc=MinInd1+MaxIntInd-1;
            

end