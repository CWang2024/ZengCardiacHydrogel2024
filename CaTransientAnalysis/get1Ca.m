function [] = get1Ca(path, filename, Output_folder_name)
%%

%path= '/Users/clarkwang/Documents/MATLAB/AP&Ca analysis/BatchCaAnalysis/videoCa analysis/para-10x.mp4';
%filename='para-10x.mp4';
%Output_folder_name='ZQJCa_Output';
baseName=filename(1:end-4);
cd(Output_folder_name);
%mkdir(baseName);
%cd(baseName);
vidObj = VideoReader(path);

% Preallocate for averaging frames
firstFrame = im2gray(readFrame(vidObj));  % Store the first frame for display
avgFrame = double(firstFrame);            % Initialize with the first frame (as double for averaging)

% Process each frame and compute average
while hasFrame(vidObj)
    frame = readFrame(vidObj);
    grayFrame = im2gray(frame);           % Convert to grayscale
    avgFrame = avgFrame + double(grayFrame);  % Sum all frames
end

% Calculate the average image
avgFrame = avgFrame / vidObj.NumFrames;

% Contrast enhancement
enhancedFrame = imadjust(uint8(avgFrame));  % Convert back to uint8 and enhance contrast
screenSize = get(0, 'ScreenSize');
squareSize = min(screenSize(3), screenSize(4));

% Adjust the image and display it in a figure

figure('Position', [0, 0, squareSize, squareSize]);
% Create subplots


% Left: Display original frame 1
%subplot(1, 2, 1);
%imshow(firstFrame);
%title('Original Frame 1');

% Right: Display contrast-enhanced averaged frame
%subplot(1, 2, 2);
imshow(enhancedFrame,'InitialMagnification', 'fit');
title('Enhanced Merged & Averaged Image');
%IntData=zeros(vidObj.NumFrames, 1);
%%

P = {};
disp("Draw ROIs, press Enter key and triple click any points on the image to end.");
%disp("The last ROI should be a dark background, not a target cell.")




%while true

    h = drawpolyline(); 
    %key = get(gcf, 'CurrentCharacter');
    wait(h);
    %{
    if key == 13
        key=0;
        break;
    elseif key == 32
        key=0;
        continue;
    end
    %}
    
    P{end+1} = h.Position;

    roiMask = poly2mask(P{end}(:, 1), P{end}(:, 2), size(frame, 1), size(frame, 2));
    stats = regionprops(roiMask, 'Centroid');
    centroid = stats.Centroid;

        % Add label at the centroid
    text(centroid(1), centroid(2), 'Ca ROI', 'Color', 'red', 'FontSize', 12);
    % Wait for key press to check if 'q' is pressed

%end

saveas(gcf, strcat(baseName, '_ROI selection.tif'));


numROI=length(P);
if numROI==0
    disp("Skip current video ... initializing next video ... ...")
    IndexOut=index;
    clear figure;
    clf;
    close all;
    return;
end


disp("ROI selection complete, running analysis...")

% Save the figure with ROI labels as PNG
%output_fig_filename = fullfile(Output_dir, [Input_file_name(1:end-4), '_ROI_labels.png']);
%output_xlsx_filename = fullfile(Output_dir, [Input_file_name(1:end-4), '.xlsx']);
%saveas(gcf, output_fig_filename);

%% 
Masks=cell(1,numROI);
for i=1:numROI+1
    if i==numROI+1
        Masks{i} = uint8(ones(size(frame, 1), size(frame, 2)));
    else
        Masks{i} = uint8(poly2mask(P{i}(:,1), P{i}(:,2), size(frame, 1), size(frame, 2)));
    end
end
vidObj = VideoReader(path);
NumFrames = vidObj.NumFrames;
IntDataDrive=zeros(NumFrames, numROI);



for i=1:NumFrames
    frame=readFrame(vidObj);
    for j=1:numROI+1
        temp=frame.*Masks{j};
        temp=temp(temp~=0);
        IntDataDrive(i,j)= mean(temp,"all");
    end
end

%%
%{
vidObj = VideoReader(path);
for i = 1:vidObj.NumFrames
    frame=readFrame(vidObj);
    frame=im2gray(frame);
    IntData(i)=mean(frame, "All");
end


L=1;
MeanDeltaInt = zeros(1, L);
MeanRatioInt=zeros(1,L);
MeanTTP10 = zeros(1, L);
MeanTTP50 = zeros(1, L);
MeanTTP90 = zeros(1, L);
MeanTTD10 = zeros(1, L);
MeanTTD50 = zeros(1, L);
MeanTTD90 = zeros(1, L);
MeanCaRelease = zeros(1,L);
RMeanCaRelease = zeros(1,L);
MeanCaIntake = zeros(1,L);
RMeanCaIntake = zeros(1,L);
Peakfrequency = zeros(1,L);
TissueIndex=strings(1,L);
MeanFilterError = zeros(1, L);
Beat2BeatInterval = zeros(1, L);

    for i = 1:L
    
        Index = num2str(i);
        TissueIndex(i)=filename(1:end-4);
        OriginalIntData = IntData;
        
      
        [p,S,mu] = polyfit([1:length(IntData)]', IntData, 2); 
        IntDataBaseline = polyval(p, [1:length(IntData)]', [], mu);
        IntData = IntData - IntDataBaseline + mean(IntData);
        IntData = lowpass(IntData,floor(fps/2)-2,fps);
        
    
        %Filter Data
        fps = vidObj.FrameRate;
        Times=(1/fps)*[1:1:vidObj.NumFrames]';
        
        OriginalTimes=Times;
        
        plot(OriginalTimes, OriginalIntData,'LineWidth', 2 );
        title(strcat(filename(1:end-4), "Raw"));
        ylabel("MFI (a.u.)");
        xlabel("Time (s)");
        plotFilename = strcat(filename(1:end-4), '_Raw.png');
        plotPath = fullfile(Output_folder_name, plotFilename);
        saveas(gcf, plotPath);
        close(gcf);
        

        
        FilterData = lowpass(IntData,floor(fps/2)-2,fps);%F needs to be lower than at least half the frame rate of video to be within Nyquist frequency
    
        IntData = FilterData(40:(end-40));
        Times = Times(40:(end-40));

        %IntData=IntData(15:(end - 15));
        
        
       % total recording, should extend the recording time longer 20-30s
       
        subplot(L,1,i);
        plot(Times,IntData, 'LineWidth', 2);
        
        title(strcat(filename(1:end-4), "Corrected"));
        ylabel("MFI (a.u.)");
        xlabel("Time (s)");
        hold on;
        %Find peaks
        MinProm = 0.5.*(max(IntData) - min(IntData));%0.1 set by joe
        [PeakValue, PeakInd] = findpeaks(IntData,'MinPeakProminence', MinProm);
        PeakTime = Times(PeakInd);
       
        %Determine minimum values for each twitch
        MinValues = zeros(length(PeakInd)-1,1);
        MinInd = zeros(length(PeakInd)-1,1);
        % MinInd = zeros(length(PeakInd)-2,1) Wrong indexing
       
        for t = 2:length(PeakValue)
            TraceInd = PeakInd(t-1):PeakInd(t);
            Trace = IntData(TraceInd);
            MinValues(t-1) = min(Trace);
            MinInd(t-1) = find(Trace == MinValues(t-1),1,'last') + PeakInd(t-1) - 1;
        end
       
        MinTime=Times(MinInd);
        scatter(PeakTime, PeakValue);
        hold on;
        scatter(MinTime, MinValues);
        hold on;
    
        DeltaInt = zeros(1,length(PeakValue)-2);
        RatioInt=zeros(1,length(PeakValue)-2);
        TTP10 = zeros(1,length(PeakValue)-2);
        TTP50 = zeros(1,length(PeakValue)-2);
        TTP90 = zeros(1,length(PeakValue)-2);
        TTD10 = zeros(1,length(PeakValue)-2);
        TTD50 = zeros(1,length(PeakValue)-2);
        TTD90 = zeros(1,length(PeakValue)-2);
        CaRelease = zeros(1,length(PeakValue)-2);
        RCaRelease = zeros(1,length(PeakValue)-2);
        CaIntake = zeros(1,length(PeakValue)-2);
        RCaIntake = zeros(1,length(PeakValue)-2);
        FilterError = zeros(1,length(PeakValue)-2);
        
        
        RMMSD=0;
        for t = 2:(length(PeakValue)-1)
    
            TraceInd = MinInd(t-1):MinInd(t);
            TraceTimes = Times(TraceInd);
            Trace = IntData(TraceInd);
            OriginalTrace=OriginalIntData(TraceInd);
    
            MaxInt = max(Trace);
            MaxIntInd = find(Trace == MaxInt,1,'first');
            
            UpProm = MaxInt - Trace(1);
            DownProm = MaxInt - Trace(end);
            StartInd=1;
            for j=1:MaxIntInd-2
                if (Trace(j+2)-Trace(j))>0.15*UpProm && Trace(j+1)-Trace(j)>0 && Trace(j+2)-Trace(j+1)>0
                   StartInd=j;
                   break;
                end
            end

            EndInd=length(Trace);
            Down80Ind = find(Trace >= Trace(end) + 0.2.*DownProm,1,'last');
            for n=Down80Ind+1:length(Trace)-4
                if (Trace(n+4)-Trace(n)>0.3*(Trace(n+1)-Trace(n)))&&(Trace(n+1)-Trace(n)<0&&Trace(n+5)-Trace(n)>0.5*(Trace(n+1)-Trace(n)))
                    EndInd=n;
                    break;
                end
            end


    
            %MinInt = MinValues(t-1);
            %MinIntInd = find(Trace == MinInt,1,'first');
            TraceUp=Trace(StartInd:MaxIntInd);
            TraceDown=Trace(MaxIntInd:EndInd);
            TimeUp=TraceTimes(StartInd:MaxIntInd);
            TimeDown=TraceTimes(MaxIntInd:EndInd);
    
            UpProm = MaxInt - Trace(StartInd);
            DownProm = MaxInt - Trace(EndInd);
    
            %Up10Ind = find(Trace >= Trace(1) + 0.1.*UpProm,1,'first');
            %Up50Ind = find(Trace >= Trace(1) + 0.5.*UpProm,1,'first');
            %Up90Ind = find(Trace >= Trace(1) + 0.9.*UpProm,1,'first');
            %Down10Ind = find(Trace >= Trace(end) + 0.9.*DownProm,1,'last');
            %Down50Ind = find(Trace >= Trace(end) + 0.5.*DownProm,1,'last');
            %Down90Ind = find(Trace >= Trace(end) + 0.1.*DownProm,1,'last');

            F0=0.5*(TraceUp(1)+TraceDown(end));
            DeltaInt(t-1) = MaxInt - F0;
            RatioInt(t-1)=DeltaInt(t-1)/F0;
            [OutTime10, OutTrace10, OutDouble10] = findUp(TraceUp, TimeUp, 0.1, UpProm, fps);
            [OutTime50, OutTrace50, OutDouble50] = findUp(TraceUp, TimeUp, 0.5, UpProm, fps);
            [OutTime90, OutTrace90, OutDouble90] = findUp(TraceUp, TimeUp, 0.9, UpProm, fps);
            [dOutTime10, dOutTrace10, dOutDouble10] = findDown(TraceDown, TimeDown, 0.1, DownProm, fps);
            [dOutTime50, dOutTrace50, dOutDouble50] = findDown(TraceDown, TimeDown, 0.5, DownProm, fps);
            [dOutTime90, dOutTrace90, dOutDouble90] = findDown(TraceDown, TimeDown, 0.9, DownProm, fps);
            %function [OutTime, OutTrace, OutDouble] = findUp(Trace, Times, ratio, UpProm, fps)
            %TTP10(t-1) = TraceTimes(MaxIntInd) - TraceTimes(Up10Ind);
            %TTP50(t-1) = TraceTimes(MaxIntInd) - TraceTimes(Up50Ind);
            %TTP90(t-1) = TraceTimes(MaxIntInd) - TraceTimes(Up90Ind);
            MaxTime=TraceTimes(MaxIntInd);
            TTP10(t-1) = MaxTime - OutTime10;
            TTP50(t-1) = MaxTime - OutTime50;
            TTP90(t-1) = MaxTime - OutTime90;            
            TTD10(t-1) = dOutTime10 - MaxTime;
            TTD50(t-1) = dOutTime50 - MaxTime;
            TTD90(t-1) = dOutTime90 - MaxTime;
            CaRelease(t-1) = 0.9*DeltaInt(t-1)/TTP10(t-1);
            RCaRelease(t-1) = 0.9*RatioInt(t-1)/TTP10(t-1);
            CaIntake(t-1) = 0.9*DeltaInt(t-1)/TTD90(t-1);
            RCaIntake(t-1) = 0.9*RatioInt(t-1)/TTD90(t-1);
  %          FilterError(t-1) =mean(abs(Trace(Up10Ind:Down90Ind)- OriginalTrace(Up10Ind:Down90Ind)));
            FilterError(t-1) =mean(abs(Trace(:)- OriginalTrace(:)));
            
            
    
%{            
            scatter(TraceTimes(Up10Ind),Trace(Up10Ind),  "r" );
            hold on;
            scatter(TraceTimes(Up50Ind),Trace(Up50Ind),  "g"   );
            hold on;
            scatter(TraceTimes(Up90Ind),Trace(Up90Ind),  "b"   );
            hold on;
%}
            scatter(OutTime10, OutTrace10,  "r" );
            hold on;
            scatter(OutTime50, OutTrace50,  "g"   );
            hold on;
            scatter(OutTime90, OutTrace90,  "b"   );
            hold on;
            scatter(dOutTime10,dOutTrace10, "r"  );
            hold on;
            scatter(dOutTime50,dOutTrace50, "g"  );
            hold on;
            scatter(dOutTime90,dOutTrace90, "b"   );
            hold on;
            scatter(TimeUp(1),TraceUp(1), "k", 'filled'  );
            hold on;
            scatter(TimeDown(end),TraceDown(end), "r", 'filled'  );
    
        end
    
        PT2=PeakTime(2:end);
        PT1=PeakTime(1:end-1);
        BInt=PT2-PT1;
        
        RMSSD = sqrt(mean((BInt-mean(BInt)).^2));

    
    
        MeanDeltaInt(i) = mean(DeltaInt);
        MeanRatioInt(i) = mean(RatioInt);
        MeanTTP10(i) = mean(TTP10);
        MeanTTP50(i) = mean(TTP50);
        MeanTTP90(i) = mean(TTP90);
        MeanTTD10(i) = mean(TTD10);
        MeanTTD50(i) = mean(TTD50);
        MeanTTD90(i) = mean(TTD90);
        MeanCaRelease(i) = mean(CaRelease);
        RMeanCaRelease(i) = mean(RCaRelease);
        MeanCaIntake(i) = mean(CaIntake);
        RMeanCaIntake(i) = mean(RCaIntake);
        Peakfrequency(i)= (length(MinInd)-1)/(Times(MinInd(end))-Times(MinInd(1)))*60;
        MeanFilterError(i)=mean(FilterError);
        Beat2BeatInterval(i) = mean(BInt);
        HRV(i) = RMSSD;
    
    end
C1 = ["MeanDeltaInt (DeltaF (a.u.))"; "MeanRatioInt (DeltaF/F0)"; "MeanTTP10 (s)"; "MeanTTP50 (s)"; "MeanTTP90 (s)";...
                               "MeanTTD10 (s)"; "MeanTTD50 (s)"; "MeanTTD90 (s)";"MeanCaRelease (a.u./s)";"RMeanCaRelease ((DeltaF/F0)/s)";"MeanCaIntake (a.u./s)";"RMeanCaIntake ((DeltaF/F0)/s)"; "Beat rate(BPM)"; "Beat to Beat Interval(s)"; "Heart Rate Varibility (s)"; "Filter Error (a.u.)"];
AllParameters = [MeanDeltaInt; MeanRatioInt; MeanTTP10; MeanTTP50; MeanTTP90;...
                               MeanTTD10; MeanTTD50; MeanTTD90;MeanCaRelease;RMeanCaRelease;MeanCaIntake;RMeanCaIntake; Peakfrequency; Beat2BeatInterval; HRV; MeanFilterError];
%AllParameters = [C1, AllParameters];
%}

%%
allROInum = numROI+1;
numPeak = zeros(1, allROInum);
MeanDeltaInt= zeros(1, allROInum); 
MeanRatioInt = zeros(1, allROInum); 
MeanTTP10 = zeros(1, allROInum); 
MeanTTP50 = zeros(1, allROInum); 
MeanTTP90 = zeros(1, allROInum);
MeanTP = zeros(1, allROInum);
MeanTTD10 = zeros(1, allROInum); 
MeanTTD50 = zeros(1, allROInum); 
MeanTTD90 = zeros(1, allROInum);
MeanTD  = zeros(1, allROInum); 
MeancaT10 = zeros(1, allROInum); 
MeancaT50 = zeros(1, allROInum); 
MeancaT90 = zeros(1, allROInum);
MeancaT = zeros(1, allROInum);

MeanFTTP10 = zeros(1, allROInum); 
MeanFTTP50 = zeros(1, allROInum); 
MeanFTTP90 = zeros(1, allROInum);
MeanFTP = zeros(1, allROInum);
MeanFTTD10 = zeros(1, allROInum); 
MeanFTTD50 = zeros(1, allROInum); 
MeanFTTD90 = zeros(1, allROInum);
MeanFTD  = zeros(1, allROInum); 
MeanFcaT10 = zeros(1, allROInum); 
MeanFcaT50 = zeros(1, allROInum); 
MeanFcaT90 = zeros(1, allROInum);
MeanFcaT = zeros(1, allROInum);

MeanBTTP10 = zeros(1, allROInum); 
MeanBTTP50 = zeros(1, allROInum); 
MeanBTTP90 = zeros(1, allROInum);
MeanBTP = zeros(1, allROInum);
MeanBTTD10 = zeros(1, allROInum); 
MeanBTTD50 = zeros(1, allROInum); 
MeanBTTD90 = zeros(1, allROInum);
MeanBTD  = zeros(1, allROInum); 
MeanBcaT10 = zeros(1, allROInum); 
MeanBcaT50 = zeros(1, allROInum); 
MeanBcaT90 = zeros(1, allROInum);
MeanBcaT = zeros(1, allROInum);

MeanTau = zeros(1, allROInum); 
MeanRsqr = zeros(1, allROInum); 
MeanCaRelease = zeros(1, allROInum);
RMeanCaRelease = zeros(1, allROInum);
MeanCaIntake = zeros(1, allROInum);
RMeanCaIntake = zeros(1, allROInum); 
BRM = zeros(1, allROInum); 
Beat2BeatInterval = zeros(1, allROInum); 
HRV = zeros(1, allROInum); 
MeanFilterError = zeros(1, allROInum);
stdF = zeros(1, allROInum); 
stdFF0 = zeros(1, allROInum); 
stdTTP10 = zeros(1, allROInum); 
stdTTP50 = zeros(1, allROInum); 
stdTTP90 = zeros(1, allROInum); 
stdTP = zeros(1, allROInum);
stdTTD10 = zeros(1, allROInum); 
stdTTD50 = zeros(1, allROInum); 
stdTTD90 = zeros(1, allROInum);  
stdTD = zeros(1, allROInum); 
stdCaT10 = zeros(1, allROInum);
stdCaT50 = zeros(1, allROInum);
stdCaT90 = zeros(1, allROInum);
stdTau = zeros(1, allROInum);
stdupV = zeros(1, allROInum); 
RstdupV = zeros(1, allROInum);
stddownV = zeros(1, allROInum);
RstddownV = zeros(1, allROInum);
FilterDataDrive = zeros(size(IntDataDrive, 1)-79, size(IntDataDrive, 2));
dFF0DataDrive = zeros(size(IntDataDrive, 1)-79, size(IntDataDrive, 2));
dFF0RawDataDrive = zeros(size(IntDataDrive, 1)-79, size(IntDataDrive, 2));
fps = vidObj.FrameRate;

for i = 1:allROInum
        IntData = IntDataDrive(:, i);
        OriginalIntData = IntData;
        
        
        %{
        [p,S,mu] = polyfit([1:length(IntData)]', IntData, 2); 
        IntDataBaseline = polyval(p, [1:length(IntData)]', [], mu);
        IntData = IntData - IntDataBaseline + mean(IntData);
        %}
        
        %IntData=smoothTrace(IntData, 11, 2);
        %{
        [p,~,mu] = polyfit((1:length(IntData))', IntData, 2); 
        IntDataBaseline = polyval(p, (1:length(IntData))', [], mu);
        IntData = IntData - IntDataBaseline + mean(IntData);
%}
        IntData = lowpass(IntData,10,fps); %floor(fps/2)-2
        %Filter Data
        
        % Suitable only for large data set with multiple peaks!!!!!!
        Times=(1/fps)*(1:1:vidObj.NumFrames)';
        
        
        OriginalTimes=Times;
        subplot(3,1,1);
        plot(Times(40:(end-40)), OriginalIntData(40:(end-40)),'LineWidth', 2 );
%        title(strcat(filename(1:end-4), "Raw"));
        ylabel("MFI (a.u.)");
        xlabel("Time (s)");
%        plotFilename = strcat(filename(1:end-4), '_Raw.png');
%        plotPath = fullfile(Output_folder_name, plotFilename);
%        saveas(gcf, plotPath);
        %close(gcf);
        
        
        
        %FilterData = lowpass(IntData,5,fps);%F needs to be lower than at least half the frame rate of video to be within Nyquist frequency
        %length(IntData)
        %length(FilteredData)
        IntData = IntData(40:(end-40));
        Times = Times(40:(end-40));
        FilterDataDrive(:,i)=IntData;
        
        %IntData=IntData(15:(end - 15));
        
        
       % total recording, should extend the recording time longer 20-30s
       
%        subplot(L,1,i);
        %figure();
        subplot(3,1,2);
        plot(Times,IntData, 'LineWidth', 2);
        
%        title(strcat(filename(1:end-4), "Corrected"));
        ylabel("MFI (a.u.)");
        xlabel("Time (s)");
        hold on;
        %Find peaks
        MinProm = 0.3.*(max(IntData) - min(IntData));%0.1 set by joe
        [PeakValue, PeakInd] = findpeaks(IntData,'MinPeakProminence', MinProm);
        PeakTime = Times(PeakInd);
       
        %Determine minimum values for each twitch
        MinValues = zeros(length(PeakInd)-1,1);
        MinInd = zeros(length(PeakInd)-1,1);
        % MinInd = zeros(length(PeakInd)-2,1) Wrong indexing
       
        for t = 2:length(PeakValue)
            TraceInd = PeakInd(t-1):PeakInd(t);
            Trace = IntData(TraceInd);
            MinValues(t-1) = min(Trace);
            MinInd(t-1) = find(Trace == MinValues(t-1),1,'last') + PeakInd(t-1) - 1;
        end
        
        
        
        
        StartTraceInd=1:PeakInd(1);
        StartTrace=IntData(StartTraceInd);
        MinInt=min(StartTrace);
        MinIntInd=find(StartTrace == MinInt,1,'first');
        UpProm=StartTrace(end)-MinInt;
        if length(StartTrace)-MinIntInd>=5 
            DiffStartTrace=diff(StartTrace);
            %UpProm=StartTrace(end)-MinInt;
            
            for n=MinIntInd:PeakInd(1)-5
                    %if (Trace(j+4)-Trace(j))>0.15*UpProm && Trace(j+1)-Trace(j)>0 && Trace(j+2)-Trace(j+1)>0
                    if all(DiffStartTrace(n:n+4)>0)&& StartTrace(n+4)-StartTrace(n)>0.2*UpProm
                       if MinIntInd ==1 && n==1
                           break;
                       end
                       MinInd=[MinIntInd; MinInd];
                       MinValues=[MinInt; MinValues];
                       %disp("B1");
                       break;
                    end
            end
        

        elseif MinInt <=IntData(MinInd(1))
            MinInd=[MinIntInd; MinInd];
            MinValues=[MinInt; MinValues];
        end
        
        EndTraceInd=PeakInd(end):length(IntData);
        EndTrace=IntData(EndTraceInd);
        MinInt=min(EndTrace);
        MinIntInd=find(EndTrace == MinInt,1,'first');
        %MinInd
        %PeakInd
        
        if MinIntInd~=length(EndTrace)&&length(EndTrace)-MinIntInd>=5
            DiffEndTrace=diff(EndTrace);
            %UpProm=EndTrace(1)-MinInt;
            for n=MinIntInd:length(EndTrace)-5
                    %if (Trace(j+4)-Trace(j))>0.15*UpProm && Trace(j+1)-Trace(j)>0 && Trace(j+2)-Trace(j+1)>0
                    if all(DiffEndTrace(n:n+4)>0)&& EndTrace(n+4)-EndTrace(n)>0.2*UpProm
                       MinInd=[MinInd; MinIntInd+PeakInd(end)-1];
                       MinValues=[MinValues; MinInt];
                       sInd = n;
                       %disp("AAAAAA");
                       break;
                    end
            end
        elseif MinInt <= IntData(MinInd(end))
            MinInd=[MinInd; MinIntInd+PeakInd(end)-1];
            MinValues=[MinValues; MinInt];
        end

        %MinInd

        %for t = 2:(length(PeakValue)-1)
        %%
        %
        waveTrace = struct('Trace', [], 'TraceUp', [], 'TraceDown', [], 'Time', [], ...
                   'TimeUp', [], 'TimeDown', [], 'Startloc', [], 'Endloc', [], 'Maxloc', []);

        MinTime=Times(MinInd);
        scatter(PeakTime, PeakValue);
        hold on;
        scatter(MinTime, MinValues);
        hold on;
        
        L=length(MinInd)-1;
        DeltaInt = zeros(1,L);
        RatioInt=zeros(1,L);
        TTP10 = zeros(1,L);
        TTP50 = zeros(1,L);
        TTP90 = zeros(1,L);
        TTD10 = zeros(1,L);
        TTD50 = zeros(1,L);
        TTD90 = zeros(1,L);
        TP = zeros(1, L);
        TD = zeros(1, L);
        tau=zeros(1, L);
        Rsqr=zeros(1, L );
        CaRelease = zeros(1,L);
        RCaRelease = zeros(1,L);
        CaIntake = zeros(1,L);
        RCaIntake = zeros(1,L);
        FilterError = zeros(1,L);
        
        RMMSD=0;
        waveBaseInd = (1:MinInd(1)-1);
        for t = 1:L
            waveTrace(t) = findWaveTrace(MinInd(t), MinInd(t+1), IntData, Times);
            waveBaseInd = [waveBaseInd , (waveTrace(t).Startloc+1: waveTrace(t).Endloc-1)];
        end
        waveBaseInd = [waveBaseInd, (MinInd(end)+1:length(IntData))];
        mask=true(1,length(IntData));
        mask(waveBaseInd)=0;
        %mask=uint8(mask);
        baselineData=IntData(mask);
        F0=mean(baselineData, "all");
        baselineTime=Times(mask);
        %{
        p=polyfit(baselineTime,baselineData, 2);
        baselinefit=polyval(p, Times);
        figure();
        %}
        
        dFF0Data=(IntData-F0)/F0;
        dFF0DataDrive(:,i)=dFF0Data;
        dFF0RawData=(OriginalIntData(40:end-40)/F0);
        dFF0RawDataDrive(:,i)=dFF0RawData;
        dFF0waveTrace=raw2dFF0(waveTrace, dFF0Data);
        for t = 1:L
            %MinInd(t-1)
            %{
            TraceInd = MinInd(t-1):MinInd(t);
            TraceTimes = Times(TraceInd);
            Trace = IntData(TraceInd);
            OriginalTrace=OriginalIntData(TraceInd);
    
            MaxInt = max(Trace);
            MaxIntInd = find(Trace == MaxInt,1,'first');
            
            UpProm = MaxInt - Trace(1);
            DownProm = MaxInt - Trace(end);
            StartInd=1;
            DiffTrace=diff(Trace);
            for j=1:MaxIntInd-26
                %if (Trace(j+4)-Trace(j))>0.15*UpProm && Trace(j+1)-Trace(j)>0 && Trace(j+2)-Trace(j+1)>0
                if all(DiffTrace(j:j+25)>0)&& Trace(j+25)-Trace(j)>0.2*UpProm
                   StartInd=j;
                   break;
                end
            end
            

            EndInd=length(Trace);
            %{
            Down80Ind = find(Trace >= Trace(end) + 0.2.*DownProm,1,'last');
            for n=Down80Ind+1:length(Trace)-4
                if (Trace(n+4)-Trace(n)>0.3*(Trace(n+1)-Trace(n)))&&(Trace(n+1)-Trace(n)<0&&Trace(n+5)-Trace(n)>0.5*(Trace(n+1)-Trace(n)))
                    EndInd=n;
                    break;
                end
            end    
            %}
           
            
                 
            MaxInd=MaxIntInd-StartInd+1;
            %MinInt = MinValues(t-1);
            %MinIntInd = find(Trace == MinInt,1,'first');
            TraceUp=Trace(StartInd:MaxIntInd);
            TraceDown=Trace(MaxIntInd:EndInd);
            TimeUp=TraceTimes(StartInd:MaxIntInd);
            TimeDown=TraceTimes(MaxIntInd:EndInd);
            %}
      

    
            %UpProm = MaxInt - Trace(waveTrace(t).Startloc);
            %DownProm = MaxInt - Trace(waveTrace(t).Endloc);  
            UpProm = IntData(waveTrace(t).Maxloc) - IntData(waveTrace(t).Startloc);
            DownProm = IntData(waveTrace(t).Maxloc) - IntData(waveTrace(t).Endloc);
            %DownProm = MaxInt - Trace(waveTrace(t).Endloc);    
            %Up10Ind = find(Trace >= Trace(1) + 0.1.*UpProm,1,'first');
            %Up50Ind = find(Trace >= Trace(1) + 0.5.*UpProm,1,'first');
            %Up90Ind = find(Trace >= Trace(1) + 0.9.*UpProm,1,'first');
            %Down10Ind = find(Trace >= Trace(end) + 0.9.*DownProm,1,'last');
            %Down50Ind = find(Trace >= Trace(end) + 0.5.*DownProm,1,'last');
            %Down90Ind = find(Trace >= Trace(end) + 0.1.*DownProm,1,'last');
            p=polyfit([waveTrace(t).Startloc, waveTrace(t).Endloc], [IntData(waveTrace(t).Startloc), IntData(waveTrace(t).Endloc)], 1);
            F0=polyval(p, waveTrace(t).Maxloc);
            %F0=0.5*(TraceUp(1)+TraceDown(end));
            DeltaInt(t) = IntData(waveTrace(t).Maxloc) - F0;
            RatioInt(t)=DeltaInt(t)/F0;
            TraceUp=waveTrace(t).TraceUp;
            TraceDown=waveTrace(t).TraceDown;
            TimeUp=waveTrace(t).TimeUp;
            TimeDown=waveTrace(t).TimeDown;
            %{
            [OutTime10, OutTrace10, OutDouble10] = findUp(TraceUp, TimeUp, 0.1, UpProm, fps);
            [OutTime50, OutTrace50, OutDouble50] = findUp(TraceUp, TimeUp, 0.5, UpProm, fps);
            [OutTime90, OutTrace90, OutDouble90] = findUp(TraceUp, TimeUp, 0.9, UpProm, fps);
            [dOutTime10, dOutTrace10, dOutDouble10] = findDown(TraceDown, TimeDown, 0.1, DownProm, fps);
            [dOutTime50, dOutTrace50, dOutDouble50] = findDown(TraceDown, TimeDown, 0.5, DownProm, fps);
            [dOutTime90, dOutTrace90, dOutDouble90] = findDown(TraceDown, TimeDown, 0.9, DownProm, fps);
            %}
            
            [OutTime10, OutTrace10] = FastintersectY0(0.1*UpProm+TraceUp(1), TimeUp, TraceUp);
            [OutTime50, OutTrace50] = FastintersectY0(0.5*UpProm+TraceUp(1), TimeUp, TraceUp);
            [OutTime90, OutTrace90] = FastintersectY0(0.9*UpProm+TraceUp(1), TimeUp, TraceUp);
          
            [dOutTime10, dOutTrace10] = FastintersectY0(0.9*DownProm+TraceDown(end), TimeDown, TraceDown);
            [dOutTime50, dOutTrace50] = FastintersectY0(0.5*DownProm+TraceDown(end), TimeDown, TraceDown);
            [dOutTime90, dOutTrace90] = FastintersectY0(0.1*DownProm+TraceDown(end), TimeDown, TraceDown);
            %function [OutTime, OutTrace, OutDouble] = findUp(Trace, Times, ratio, UpProm, fps)
            %TTP10(t) = TraceTimes(MaxIntInd) - TraceTimes(Up10Ind);
            %TTP50(t) = TraceTimes(MaxIntInd) - TraceTimes(Up50Ind);
            %TTP90(t) = TraceTimes(MaxIntInd) - TraceTimes(Up90Ind);
            MaxTime=Times(waveTrace(t).Maxloc);
            TTP10(t) = MaxTime - OutTime10;
            TTP50(t) = MaxTime - OutTime50;
            TTP90(t) = MaxTime - OutTime90;            
            TTD10(t) = dOutTime10 - MaxTime;
            TTD50(t) = dOutTime50 - MaxTime;
            TTD90(t) = dOutTime90 - MaxTime;
            TP(t)= TimeUp(end)-TimeUp(1);
            TD(t)= TimeDown(end)-TimeDown(1);
            CaRelease(t) = 0.9*DeltaInt(t)/TTP10(t);
            RCaRelease(t) = 0.9*RatioInt(t)/TTP10(t);
            CaIntake(t) = 0.9*DeltaInt(t)/TTD90(t);
            RCaIntake(t) = 0.9*RatioInt(t)/TTD90(t);
            [tau(t), Rsqr(t)] = getTimeConstant(TimeDown, TraceDown);
            if Rsqr(t)<0.9
                tau(t) = NaN;
                Rsqr(t) = NaN;
            end
  %          FilterError(t) =mean(abs(Trace(Up10Ind:Down90Ind)- OriginalTrace(Up10Ind:Down90Ind)));
            FilterError(t) =mean(abs(OriginalIntData(waveTrace(t).Startloc:waveTrace(t).Endloc)-waveTrace(t).Trace));

            scatter(OutTime10, OutTrace10,  "r" );
            hold on;
            scatter(OutTime50, OutTrace50,  "g"   );
            hold on;
            scatter(OutTime90, OutTrace90,  "b"   );
            hold on;
            scatter(dOutTime10,dOutTrace10, "r"  );
            hold on;
            scatter(dOutTime50,dOutTrace50, "g"  );
            hold on;
            scatter(dOutTime90,dOutTrace90, "b"   );
            hold on;
            scatter(TimeUp(1),TraceUp(1), "k", 'filled'  );
            hold on;
            scatter(TimeDown(end),TraceDown(end), "r", 'filled'  );
            %hold on;
        
        end
        caT10=TTP90+TTD10;
        caT50=TTP50+TTD50;
        caT90=TTP10+TTD90;
        caT = TP + TD;

        PT2=PeakTime(2:end);
        PT1=PeakTime(1:end-1);
        BInt=PT2-PT1;
        
        RMSSD = sqrt(mean((BInt-mean(BInt)).^2));
        
    
    
        MeanDeltaInt(i) = mean(DeltaInt);
        MeanRatioInt(i) = mean(RatioInt);
        MeanTTP10(i) = mean(TTP10);
        MeanTTP50(i) = mean(TTP50);
        MeanTTP90(i) = mean(TTP90);
        MeanTTD10(i) = mean(TTD10);
        MeanTTD50(i) = mean(TTD50);
        MeanTTD90(i) = mean(TTD90);
        MeancaT10(i) = mean(caT10);
        MeancaT50(i) = mean(caT50);
        MeancaT90(i) = mean(caT90);
        MeanTD(i) = mean(TD);
        MeanTP(i) = mean(TP);
        MeancaT(i) = mean(caT);

        

        MeanTau(i) =mean(tau, 'omitnan')*1000;
        MeanRsqr(i) =mean(Rsqr, 'omitnan');
        MeanCaRelease(i)  = mean(CaRelease);
        RMeanCaRelease(i)  = mean(RCaRelease);
        MeanCaIntake(i)  = mean(CaIntake);
        RMeanCaIntake(i)  = mean(RCaIntake);
        BRM(i) = L/(Times(MinInd(end))-Times(MinInd(1)))*60;
        if length(PeakInd)>= 2
            BRM(i) = 60/mean(diff(PeakTime));
        end

        MeanFilterError(i) = mean(FilterError);
        Beat2BeatInterval(i) = mean(BInt);
        HRV(i) = RMSSD;
        numPeak(i) = L;
%%      
        FR = sqrt(Beat2BeatInterval(i));
        MeanFTTP10(i) = MeanTTP10(i)/FR; 
        MeanFTTP50(i) = MeanTTP50(i)/FR; 
        MeanFTTP90(i) = MeanTTP90(i)/FR;
        MeanFTP(i) = MeanTP(i)/FR;
        MeanFTTD10(i) = MeanTTD10(i)/FR; 
        MeanFTTD50(i) = MeanTTD50(i)/FR; 
        MeanFTTD90(i) = MeanTTD90(i)/FR;
        MeanFTD(i)  = MeanTD(i)/FR; 
        MeanFcaT10(i) = MeancaT10(i)/FR; 
        MeanFcaT50(i) = MeancaT50(i)/FR; 
        MeanFcaT90(i) = MeancaT90(i)/FR;
        MeanFcaT(i) = MeancaT(i)/FR;
        
        BR = nthroot(Beat2BeatInterval(i), 3);
        MeanBTTP10(i) = MeanTTP10(i)/BR; 
        MeanBTTP50(i) = MeanTTP50(i)/BR; 
        MeanBTTP90(i) = MeanTTP90(i)/BR;
        MeanBTP(i) = MeanTP(i)/BR;
        MeanBTTD10(i) = MeanTTD10(i)/BR; 
        MeanBTTD50(i) = MeanTTD50(i)/BR; 
        MeanBTTD90(i) = MeanTTD90(i)/BR;
        MeanBTD(i)  = MeanTD(i)/BR; 
        MeanBcaT10(i) = MeancaT10(i)/BR; 
        MeanBcaT50(i) = MeancaT50(i)/BR; 
        MeanBcaT90(i) = MeancaT90(i)/BR;
        MeanBcaT(i) = MeancaT(i)/BR;

        hold off;
        subplot(3, 1, 3);
        plot(Times, dFF0Data, "r", 'LineWidth', 2);
       
%       title(strcat(filename(1:end-4), "Raw"));
        ylabel("dF/F");
        xlabel("Time (s)");
        
        
        %figure();
        %plot(Times, dFF0Data);
        if i==1

            avgTrace = calculateAverageTrace(waveTrace)';
            avgTime = (1:1:length(avgTrace))/fps;
            avgTime = avgTime';
            dFF0avgTrace = calculateAverageTrace(dFF0waveTrace)';
        end
        %plotAll(truncatedTraces);
        
        
        stdF(i) = std(DeltaInt);
        stdFF0(i) =std(RatioInt); 
        stdTTP10(i) = std(TTP10); 
        stdTTP50(i) = std(TTP50); 
        stdTTP90(i) = std(TTP90); 
        stdTP(i) = std(TP);
        stdTTD10(i)= std(TTD10); 
        stdTTD50(i) =std(TTD50); 
        stdTTD90(i) =std(TTP90);  
        stdTD(i) = std(TD); 
        stdCaT10(i) = std(caT10);
        stdCaT50(i) = std(caT50);
        stdCaT90(i) = std(caT90);
        stdTau(i) = std(tau, 'omitnan')*1000;
        stdupV(i) = std(CaRelease); 
        RstdupV(i) = std(RCaRelease);
        stddownV(i) = std(CaIntake);
        RstddownV(i) = std(RCaIntake);
        if i==numROI+1
            plotFilename = strcat('WholeImage_', filename(1:end-4), '.tif');
        else
            plotFilename = strcat('SelectedROI', '_', filename(1:end-4), '.tif');
        end
        %plotPath = fullfile(Output_folder_name, plotFilename);
        saveas(gcf, plotFilename);
        close(gcf);
end
%%   
C1 = ["Peak number analyzed"; "ΔF (a.u.)"; "ΔF/F0 "; "TTP10 (s)"; "TTP50 (s)"; "TTP90(s)"; "TP(s)";...
                               "TTD10 (s)"; "TTD50 (s)"; "TTD90 (s)";  "TD(s)";"CaT10(s)";"CaT50(s)";"CaT90(s)"; "CaT(s)"; "F_TTP10 (s)"; "F_TTP50 (s)"; "F_TTP90(s)"; "F_TP(s)";...
                               "F_TTD10 (s)"; "F_TTD50 (s)"; "F_TTD90 (s)";  "F_TD(s)";"F_CaT10(s)";"F_CaT50(s)";"F_CaT90(s)"; "F_CaT(s)";"B_TTP10 (s)"; "B_TTP50 (s)"; "B_TTP90(s)"; "B_TP(s)";...
                               "B_TTD10 (s)"; "B_TTD50 (s)"; "B_TTD90 (s)";  "B_TD(s)";"B_CaT10(s)";"B_CaT50(s)";"B_CaT90(s)"; "B_CaT(s)"; "𝜏 (ms)"; "R²";...
                               "Absolute UpSpeed (a.u./s) "; "Relative UpSpeed ((DeltaF/F0)/s)";"Absolute DownSpeed (a.u./s)";"Relative DownSpeed ((DeltaF/F0)/s)";...
                               "Beats/min"; "Beat to Beat Interval(s)"; "Beat Interval Variability(s)"; "Filter Error (a.u.)";...
                               "stdΔF"; "stdΔF/F0 "; "stdTTP10 (s)"; "stdTTP50 (s)"; "stdTTP90 (s)"; "stdTP (s)";...
                               "stdTTD10 (s)"; "stdTTD50 (s)"; "stdTTD90 (s)";  "stdTD (s)";"stdCaT10";"stdCaT50";"stdCaT90";"std Time Constant";...
                               "std AbsoluteUpSpeed (a.u./s) "; "std RelativeUpSpeed ((ΔF/F0)/s)";"std AbsoluteDownSpeed (a.u./s)";"std RelativeDownSpeed ((ΔF/F0)/s)"];
AllParameters = [numPeak; MeanDeltaInt; MeanRatioInt; MeanTTP10; MeanTTP50; MeanTTP90;MeanTP;...
                               MeanTTD10; MeanTTD50; MeanTTD90;MeanTD; MeancaT10; MeancaT50; MeancaT90; MeancaT;...
                               MeanFTTP10; MeanFTTP50; MeanFTTP90;MeanFTP;...
                               MeanFTTD10; MeanFTTD50; MeanFTTD90;MeanFTD; MeanFcaT10; MeanFcaT50; MeanFcaT90; MeanFcaT;...
                               MeanBTTP10; MeanBTTP50; MeanBTTP90;MeanBTP;...
                               MeanBTTD10; MeanBTTD50; MeanBTTD90;MeanBTD; MeanBcaT10; MeanBcaT50; MeanBcaT90; MeanBcaT;...
                               MeanTau; MeanRsqr; MeanCaRelease;RMeanCaRelease;MeanCaIntake;RMeanCaIntake; BRM; Beat2BeatInterval; HRV; MeanFilterError;...
                               stdF; stdFF0; stdTTP10; stdTTP50; stdTTP90; stdTP;...
                               stdTTD10; stdTTD50; stdTTD90;  stdTD; stdCaT10;stdCaT50;stdCaT90;stdTau;...
                               stdupV; RstdupV;stddownV;RstddownV];
%{
numPeak, MeanDeltaInt, MeanRatioInt, MeanTTP10, MeanTTP50, MeanTTP90,MeanTP,...
                               MeanTTD10, MeanTTD50, MeanTTD90,MeanTD, MeancaT10, MeancaT50, MeancaT90, MeancaT...
                               MeanFTTP10, MeanFTTP50, MeanFTTP90,MeanFTP,...
                               MeanFTTD10, MeanFTTD50, MeanFTTD90,MeanFTD, MeanFcaT10, MeanFcaT50, MeanFcaT90, MeanFcaT...
                               MeanBTTP10, MeanBTTP50, MeanBTTP90,MeanBTP,...
                               MeanBTTD10, MeanBTTD50, MeanBTTD90,MeanBTD, MeanBcaT10, MeanBcaT50, MeanBcaT90, MeanBcaT...
                               MeanTau, MeanRsqr, MeanCaRelease,RMeanCaRelease,MeanCaIntake,RMeanCaIntake, BRM, Beat2BeatInterval, HRV, MeanFilterError,...
                               stdF, stdFF0, stdTTP10, stdTTP50, stdTTP90, stdTP,...
                               stdTTD10, stdTTD50, stdTTD90,  stdTD, stdCaT10,stdCaT50,stdCaT90,stdTau,...
                               stdupV, RstdupV,stddownV,RstddownV],

        %plotFilename = strcat(filename(1:end-4), '_Corrected.tif');
        %plotPath = fullfile(Output_folder_name, plotFilename);
        %saveas(gcf, plotPath);
%} 

%%
Times=(1/fps)*(1:1:vidObj.NumFrames)';
Times2=Times(40:end-40);
        excelFilename = strcat(baseName, '.xlsx');
        excelPath=excelFilename;
        %excelPath = fullfile(Output_folder_name, excelFilename);
        headerBase = cell(1,allROInum);
     
        for i=1:allROInum
            if i==allROInum
             headerBase{i} = 'Whole Image';
            else
             headerBase{i} = strcat('ROI_', num2str(i));
            end
        end

        
        % Write AllParameters to sheet 1
        %writematrix(AllParameters, excelPath, 'Sheet', 'All Parameters');
        headers = [{'All Parameters'}, headerBase];
        CaParameters = [C1, num2cell(AllParameters)];
        ParameterData = [headers; CaParameters];
        writematrix(ParameterData, excelPath, 'Sheet', 1);
    
        
        % Write OriginalTimes and OriginalIntData to sheet 2
        %writematrix([OriginalTimes, OriginalIntData], excelPath, 'Sheet', 'Raw Trace');
        % Prepare the data with headers for "Raw Trace"
        
        headers= [{'Raw Trace _ Time'}, headerBase];
        rawTraceData = [headers; num2cell([Times, IntDataDrive])];
        
        % Write the data with headers to sheet "Raw Trace"
        writecell(rawTraceData, excelPath, 'Sheet', 'Raw dF');
        
        % Write Times and IntData to sheet 3
        
        headers= [{'Filtered Trace _ Time'}, headerBase];
        TraceData = [headers; num2cell([Times2, FilterDataDrive])];
        
        % Write the data with headers to sheet "Raw Trace"
        writecell(TraceData, excelPath, 'Sheet', 'Filtered dF');
        
        headers= [{'ΔF/F0 _ Time'}, headerBase];
        TraceData = [headers; num2cell([Times2, dFF0RawDataDrive])];
        
        % Write the data with headers to sheet "Raw Trace"
        writecell(TraceData, excelPath, 'Sheet', 'dFF0');

        headers= [{'Filtered ΔF/F0 _ Time'}, headerBase];
        dFData = [headers; num2cell([Times2, dFF0DataDrive])];
        
        % Write the data with headers to sheet "Raw Trace"
        writecell(dFData, excelPath, 'Sheet', 'Filtered dFF0');
        %writematrix([Times, IntData], excelPath, 'Sheet', 'Filtered Trace');
        headers= [{'AvgTime'}, {'AvgTrace'}];
        
        AvgData = [headers; num2cell([avgTime, avgTrace])];
        
        % Write the data with headers to sheet "Raw Trace"
        writecell(AvgData, excelPath, 'Sheet', 'AvgTrace');

        dFF0AvgData = [headers; num2cell([avgTime, dFF0avgTrace])];
        % Write the data with headers to sheet "Raw Trace"
        writecell(dFF0AvgData, excelPath, 'Sheet', 'dFF0AvgTrace');
        %end
        %cd ..
        cd ..

       
end



function avgTrace = calculateAverageTrace(waveTrace)
    % Find the shortest trace length
    maxLength = max(arrayfun(@(x) length(x.Trace), waveTrace));
    
    % Initialize a matrix to store all truncated traces
    truncatedTraces = ones(length(waveTrace), maxLength);%minlength
    
    % Loop through each waveTrace element and truncate the trace to the minimum length
    for i = 1:length(waveTrace)
        dist=length(waveTrace(i).Trace);
        truncatedTraces(i, 1:dist) = waveTrace(i).Trace(1:end)';
        if dist<maxLength
            truncatedTraces(i, dist + 1: end) = waveTrace(i).Trace(end)*truncatedTraces(i, dist + 1: end);
        end
    end
    
    % Calculate the average trace
    avgTrace = mean(truncatedTraces, 1);
end

function [] = plotAll(allTrace)
    figure();
    for i=1:(size(allTrace, 1))
        plot(allTrace(i, :));
        hold on;
    end
    hold off;
end

function dFF0waveTrace=raw2dFF0(waveTrace, dFF0Data)
    dFF0waveTrace=waveTrace;
    for i=1:length(waveTrace)
        dFF0waveTrace(i).Trace=dFF0Data(waveTrace(i).Startloc:waveTrace(i).Endloc);
        dFF0waveTrace(i).TraceUp=dFF0Data(waveTrace(i).Startloc:waveTrace(i).Maxloc);
        dFF0waveTrace(i).TraceDown=dFF0Data(waveTrace(i).Maxloc:waveTrace(i).Endloc);
    end
end

function [tau, R2] = getTimeConstant(TimeDown, TraceDown)
    % Define the exponential decay function to fit
    expDecay = @(coeffs, x) coeffs(1) * exp(coeffs(2) * x) + coeffs(3);

    % Initial guess for the parameters [a, b, c]
    % a: amplitude of the decay
    % b: decay rate (which will be used to calculate tau)
    % c: offset (baseline)
    initialGuess = [max(TraceDown), -1, min(TraceDown)];
    TemptD=TimeDown-TimeDown(1);

    % Perform the curve fitting
    options = optimset('Display', 'off'); % Suppress output
    coeffs = lsqcurvefit(expDecay, initialGuess, TemptD, TraceDown, [], [], options);

    % Extract the decay rate (b) to calculate tau
    b = coeffs(2);
    tau = -1 / b;
    
    % Generate the fitted curve using the found coefficients
    TraceFit = expDecay(coeffs, TemptD);
    SS_res = sum((TraceDown - TraceFit).^2);    % Residual sum of squares›
    SS_tot = sum((TraceDown - mean(TraceDown)).^2); % Total sum of squares
    R2 = 1 - (SS_res / SS_tot);   % R-squared formula
    
    % Plot the original decay trace and the fitted curve
    
    %plot(TimeDown, TraceFit, 'r-', 'DisplayName', sprintf('Fitted Curve: y = %.2f * exp(%.2f * x) + %.2f', coeffs(1), coeffs(2), coeffs(3)));
    plot(TimeDown, TraceFit, 'r-');
    hold on;
end


function [xout, yout] = FastintersectY0(y0, x, y)

    % INTERSECT_WITH_Y0 Finds intersections of a curve with a horizontal line y = y0
    %   [XOUT, YOUT] = INTERSECT_WITH_Y0(Y0, X, Y) computes the points of 
    %   intersection between the curve defined by X and Y, and the horizontal line y = Y0.
    %
    %   X and Y are vectors where each pair (X(i), Y(i)) defines a point on the curve.
    %   XOUT and YOUT are the x and y coordinates of the intersection points.

    % Initialize output vectors
    xout = [];
    %yout = [];

    % Ensure x and y are column vectors
    x = x(:);
    %y = y(:);
    %startInd=find(y>y0, 1, "first")-1;
    %EndInd=find(y<y0, 1, "last")+1;
    diff=y-y0;
    % Loop through each segment of the curve
    for i = 1:length(y)-1
        % Check if the curve crosses y = y0 between points i and i+1
        if (diff(i)) * (diff(i+1)) < 0
            % Linear interpolation to find the exact intersection
            x1 = x(i);
            x2 = x(i + 1);
            y1 = y(i);
            y2 = y(i + 1);
            
            % Compute the x-coordinate of the intersection
            x_intersect = x1 + (y0 - y1) * (x2 - x1) / (y2 - y1);
            %y_intersect = y0;

            % Append to output
            xout(end+1, 1) = x_intersect;
            %yout(end+1, 1) = y_intersect;
        elseif y(i) == y0
            % If a point is exactly on the line y = y0
            xout(end+1, 1) = x(i);
            %yout(end+1, 1) = y(i);
        end
    end
    
    % Check the last point if it's exactly on the line y = y0
    if y(end) == y0
        xout(end+1, 1) = x(end);
        %yout(end+1, 1) = y(end);
    end
    

    xout=mean(xout);
    yout=y0;
end

function waveTrace = findWaveTrace(MinInd1, MinInd2, IntData, Times)


% Start point is defined as the point with sharp increase, end point is
% the point that the MFI is equal to the start point MFI, if no such point
% exist, then the min point is used. 
           
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
            for j=1:MaxIntInd-5
                %if (Trace(j+4)-Trace(j))>0.15*UpProm && Trace(j+1)-Trace(j)>0 && Trace(j+2)-Trace(j+1)>0
                if all(DiffTrace(j:j+4)>0)&& Trace(j+4)-Trace(j)>0.2*UpProm
                   StartInd=j;
                   break;
                end
            end
            
            EndTrace=Trace(MaxIntInd:end);
            EndInd=find(EndTrace <= Trace(StartInd),1,'first');
            if isempty(EndInd)
                MinInt=min(EndTrace);
                EndInd = find(EndTrace == MinInt,1,'first');
            end 
            EndInd=EndInd+MaxIntInd-1;
            
            %{
            Down70Ind = find(Trace >= Trace(end) + 0.3*DownProm,1,'last');
            for n=Down70Ind+2:length(Trace)-3
                
                if (Trace(n+2)-Trace(n)>0.2*(Trace(n)-Trace(n-2)))&&(DiffTrace(n)>-0.10)&&DiffTrace(n)<DiffTrace(n+1)
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