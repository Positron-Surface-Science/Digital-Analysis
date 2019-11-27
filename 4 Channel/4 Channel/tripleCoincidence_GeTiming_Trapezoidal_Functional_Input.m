function [VoutMax,timeOfFlight,errorOccurred,noFilesFound,dataFilter,baselineAverage] = ...
    tripleCoincidence_GeTiming_Trapezoidal_Functional_Input(selectPathIn,...
    iIn,analysisTypeIn,multiStopIn,multiStopConditionsIn,TypeIn,timingTypeIn,shapingTypeIn,...
    FFTIn,coinCh,ggc,RTFLIn,RTFUpIn,TFTopIn,TFRFIn,TFTopInFast,TFRFInFast,BPFSwitch,dataFilter,...
    DiscriminationVector,ELETMatrix,ClusteringNetwork,neuron,ELETParam)

%---------------------------------------------------------------------------
% Primary coincidence function. Reads the waveforms from the Lecroy
% oscilloscope files, then sends them to paes and trapezoidalfilter to
% calculate the time of flight and gamma values.
%---------------------------------------------------------------------------

selectPath = selectPathIn;%'J:\Lecroy\Coi';%uigetdir;

numChannels = 4;

channelNumber = cell(numChannels);
pulse = cell(1,4);
crossTime = cell(1,2);
errorOccurred = false;
seeingNoFiles = 0;
noFilesFound = false;
goodCounter = 0;
activeNeuron = 0;
exact = ones(size(DiscriminationVector));
if isempty(neuron)
    neuron = 0;
end
neuron = 0;
if isstruct(ELETMatrix)
    a = fieldnames(ELETMatrix);
    ELETMatrix = getfield(ELETMatrix,a{1});
end
if isempty(ggc)
    ggc = 0;
end
ggc = 1;
timingTypeChanged = false;

% PAES analysis definitions.

if multiStopIn ~= '1'
    timeOfFlight = zeros(1,10) + 1;
    
end

VoutMax = cell(1,4);
noGood = zeros(4,10);
noGoodS = zeros(4,10);
s = 1;
baselineAverage = 0;

for n=iIn:iIn+9
    
    %evalin('base','flags = 0;');
    for c=1:numChannels
        
        if TypeIn(c) ~= 4
            %'TYPEIN'
            %TypeIn(c)
            %timingTypeIn(c)
            if TypeIn(c) == 3 && (timingTypeIn(c) == 4 || ...
                    timingTypeChanged == true)
                ANN = 1;
                timingTypeIn(c) = 0;
                timingTypeChanged = true;
                'ANN ACTIVE'
            else
                ANN = 0;
                activeNeuron = 0;
                'ANN NOT ACTIVE'
            end
            %s
            %'channel #'
            %c
            noFile = true;
            noFileCounter = 1;
            noGood(c,s) = false;
            
            %---------------------------------------------------------------------------------------
            % Calling ReadLeCroyBinaryWaveform function from waveform function file to
            % extract pulse waveforms and the corresponding information into a cell of structures.
            %---------------------------------------------------------------------------------------
            
            if (true)
                channelNumber{c} = ['C',stringconversion(c),'1'];
                
                if n < 10
                    appendage = ['0000',stringconversion(n)];
                elseif n >= 10 && n < 100
                    appendage = ['000',stringconversion(n)];
                elseif n >= 100 && n < 1000
                    appendage = ['00',stringconversion(n)];
                elseif n >= 1000 && n < 10000
                    appendage = ['0',stringconversion(n)];
                else
                    appendage = stringconversion(n);
                end
                
            end
            
            if (false)
                channelNumber{c} = ['C',stringconversion(c),'Trace'];
                
                if n < 10
                    appendage = ['0000',stringconversion(n)];
                elseif n >= 10 && n < 100
                    appendage = ['000',stringconversion(n)];
                elseif n >= 100 && n < 1000
                    appendage = ['00',stringconversion(n)];
                elseif n >= 1000 && n < 10000
                    appendage = ['0',stringconversion(n)];
                else
                    appendage = stringconversion(n);
                end
                
            end
            
            if (false)
                channelNumber{c} = ['C',stringconversion(c),'XX'];
                %appendage = num2str(i-1,'%05.f');
                
                if n < 10
                    appendage = ['000000000',stringconversion(n)];
                elseif n >= 10 && n < 100
                    appendage = ['00000000',stringconversion(n)];
                elseif n >= 100 && n < 1000
                    appendage = ['0000000',stringconversion(n)];
                elseif n >= 1000 && n < 10000
                    appendage = ['000000',stringconversion(n)];
                elseif n >= 10000 && n < 100000
                    appendage = ['00000',stringconversion(n)];
                elseif n >= 100000 && n < 1000000
                    appendage = ['0000',stringconversion(n)];
                elseif n >= 1000000 && n < 10000000
                    appendage = ['000',stringconversion(n)];
                elseif n >= 10000000 && n < 100000000
                    appendage = ['00',stringconversion(n)];
                elseif n >= 10000000 && n < 100000000
                    appendage = ['0',stringconversion(n)];
                else
                    appendage = stringconversion(n);
                end
                
            end
            
            fileName = [selectPath,'\',channelNumber{c},appendage,'.trc'];
            
            % try to read the file; if it doesn't exist, try again twenty
            % times every five seconds, then break and try to read the
            % next file if its existence never occurs.
            
            while noFile == true && noFileCounter <= 7
                
                try
                    pulse{c}(s) = waveform.ReadLeCroyBinaryWaveform(fileName);
                    noFile = false;
                    
                catch
                    %errordlg('No Lecroy files within folder.','Incorrect Folder.');
                    noFile = true;
                    noFileCounter = noFileCounter + 1;
                    pause(1);
                    
                end
                
                if noFileCounter >= 5
                    noGood(c,s) = true;
                    if multiStopIn == '0'
                        timeOfFlight(c) = NaN;
                        'nofilecounter'
                        
                    else
                        timeOfFlight = zeros(1,10) + 1;
                        
                    end
                    
                    noFilesFound = true;
                    'no files'
                    return;
                    
                end
                
            end
            
            %{
            if noFile == true || isempty(pulse{c}(s).y)
                noGood(c,i) = true;
                timeOfFlight(c) = 0;
                
                if seeingNoFiles == 9
                    noFilesFound = true;
                    
                end
                
                return;
                
            end
            %}
            
            offset = -pulse{c}(s).info.OFFSET;
            Fs = pulse{c}(s).desc.fs;
            T = (1/Fs);
            oldPulse = pulse{c}(s).y;
            
            %---------------------------------------------------------------------------------------
            % HPGe pulse analysis (trapezoidal).
            %---------------------------------------------------------------------------------------
            
            if (TypeIn(c) == 3 || TypeIn(c) == 2) && shapingTypeIn(c) ~= 8
                [TPA,~,noGoodS(c,s)] = trapezoidalFilter(pulse{c}(s),oldPulse,shapingTypeIn(c),TFRFIn(c),TFTopIn(c));
                
                VoutMax{c}(s) = TPA;
                
                if noGoodS(c,s) == false
                    goodCounter = goodCounter + 1;
                    %baselineAverage = baselineAverage + baseline;
                    
                end
                
            end
            
            %---------------------------------------------------------------------------------------
            % FFT low pass filter to remove high frequency noise (on both
            % MCP and Gamma detector).
            %---------------------------------------------------------------------------------------
            
            try
                if BPFSwitch == '1'
                    [pulse{c}(s).y,dataFilter,baseline] = bpfilter(pulse{c}(s),c,FFTIn,dataFilter);
                    baseline = mean(pulse{c}(s).y(round((offset-999E-9)/T):round((offset-500E-9)/T)));
                    pulse{c}(s).y = pulse{c}(s).y - baseline;
                    
                else
                    baseline = mean(pulse{c}(s).y(round((offset-999E-9)/T):round((offset-500E-9)/T)));
                    pulse{c}(s).y = pulse{c}(s).y - baseline;
                    
                end
            catch
                'bpf error'
            end
            
            %{
            if TypeIn(1) == 4 && TypeIn(c) == 3
                assignin('base','trace',pulse{c}(s).y(round((offset-999E-9)/T):round((offset+500E-9)/T)));
                evalin('base','NaITraces = horzcat(NaITraces,trace(2550:3000));');
            
            end
            %}
            
            %---------------------------------------------------------------------------------------
            % ANN Analysis.
            %---------------------------------------------------------------------------------------
            
            %try
            if ANN == 1 && TypeIn(c) == 3
                'ANN ON'
                p = 1;
                q = floor((1.25*p)/0.5);
                FsN = (p/q)*Fs;
                TsN = 1/FsN;
                
                baseline = mean(oldPulse(round((offset-999E-9)/T):round((offset-500E-9)/T)));
                oldPulse = oldPulse - baseline;
                
                %if evalin('base','flags;') == 1
                %assignin('base','trace',pulse{c}(s).y(18626:19063));
                %evalin('base','traces = horzcat(traces,trace);');
                %end
                %{
                % ---------------------------------
                % HPGe
                vIn = resample(oldPulse,p,q);
                
                netTrace = vIn(round(14.5E-6/TsN):round(18E-6/TsN));
                
                netTrace = netTrace(150:550);
                
                netTrace = netTrace/max(netTrace);
                % ---------------------------------
                %}
                
                % ---------------------------------
                % HPGe 9x25 Network
                netTrace = pulse{c}(s).y(18626:19063);
                
                netTrace = netTrace/max(netTrace);
                % ---------------------------------
                %assignin('base','p',pulse{c}(s).y(18626:19063)/max(pulse{c}(s).y(18626:19063)));
                %crossing(s) = evalin('base','net3(p);')
                %{
                % ---------------------------------
                % NaI
                netTrace = oldPulse(round((offset-999E-9)/T):round((offset+500E-9)/T));
                
                netTrace = netTrace(2550:3000);
                netTrace = netTrace(1:250);
                
                netTrace = netTrace/min(netTrace);
                % ---------------------------------
                %}
                %assignin('base','netTrace',netTrace);
                a = fieldnames(ClusteringNetwork);
                net = getfield(ClusteringNetwork,a{1});
                a = fieldnames(DiscriminationVector);
                exact = getfield(DiscriminationVector,a{1});
                
                out = net(netTrace);
                
                'Active Neuron'
                activeNeuron = find(out == 1)
                neuron
                
                'ELET PARAMETERS'
                ELETParam
                %best = evalin('base','best;');
                ELETMatrix(activeNeuron,1:2) = ELETParam;%best(1:2,activeNeuron);
                %best(activeNeuron,1:2)
                
                'FWHM';
                %ELETMatrix(activeNeuron,4)
                
                'exact'
                exact(activeNeuron)
                
                %ELETMatrix(activeNeuron,1:2) = best(activeNeuron,1:2);
                
                if  ggc == 1 && (isempty(activeNeuron) || ...
                        isempty(neuron) || activeNeuron ~= neuron)
                    %activeNeuron = 0;
                    
                elseif ggc == 1 && isempty(activeNeuron) == 0 && ...
                        activeNeuron == neuron
                    %ELETMatrix(activeNeuron,1:2) = ELETParam;
                    %ELETMatrix(activeNeuron,3) = 0;
                    'EQUAL'
                    
                end
                
                %{
                if activeNeuron == neuron
                    plot(netTrace)
                    doesTimesExist = evalin('base','exist(''times'',''var'');');
                    
                    if doesTimesExist == 1
                        assignin('base','neuron',neuron);
                        evalin('base','times{neuron} = times{neuron} + 1;');
                        
                    else
                        evalin('base','times = cell(zeros(1,1))');
                        
                    end
                    
                end
                %}
                
            end
            %catch
            %    end
            %---------------------------------------------------------------------------------------
            % PAES analysis.
            %---------------------------------------------------------------------------------------
            
            if analysisTypeIn ~= 3 && timingTypeIn(c) ~= 3 && (c == coinCh(1) || c == coinCh(2))
                
                %{
                try
                    T = pulse{c}(s).desc.Ts;
                    offset = -pulse{c}(s).info.OFFSET;
                    pulse{c}(s).y = pulse{c}(s).y - mean(pulse{c}(s).y(round((offset-2E-6)/T):round((offset-1E-6)/T)));
                catch
                    'baseline subtraction'
                end
                %}
                % Trapezoidal shaping for timing.
                
                if TypeIn(c) == 3 && analysisTypeIn ~= 3 && timingTypeIn(c) == 2
                    [~,vOut] = trapezoidalFilter(pulse{c}(s),1,TFTopInFast,TFRFInFast);
                    pulse{c}(s).y = vOut;
                    
                end
                
                numPeaks = 1;
                [crossTimeOut,noGood(c,s),numPeaks,~] = paes_unedited(timingTypeIn(c),TypeIn(c),c,multiStopIn,...
                    multiStopConditionsIn,numPeaks,pulse{c}(s),RTFLIn,RTFUpIn,ELETMatrix, ...
                    DiscriminationVector,neuron,activeNeuron,ANN,ggc);
                
                %try
                
                    if sum(TypeIn(:) == 1) ~= 0
                        
                        switch TypeIn(c)
                            case 1
                                crossTime{s}.mcp = crossTimeOut;
                                %'cross time mcp'
                                %crossTime{s}.mcp(1)
                                
                            case 3
                                crossTime{s}.gamma = crossTimeOut;
                                
                                %{
                                if ANN == 1 && activeNeuron ~= 0 && isempty(activeNeuron) == 0
                                    'CALCULATE'
                                    crossTime{s}.gamma
                                    %crossTime{s}.gamma = crossTime{s}.gamma - ELETMatrix(activeNeuron,3);
                                    crossTime{s}.gamma
                                end
                                %}
                                
                                %'cross time gamma'
                                %crossTime{s}.gamma(1)
                        end
                        
                    elseif c == 1
                        crossTime{s}.gamma1 = min(abs(pulse{c}(s).x));
                        
                    elseif c == 2
                        crossTime{s}.gamma2 = crossTimeOut;
                        
                    elseif c == 3
                        crossTime{s}.gamma3 = crossTimeOut;
                        
                    elseif c == 4
                        crossTime{s}.gamma4 = crossTimeOut;
                        
                    end
               %{     
                catch
                    errordlg('Check channel timing type.','Timing Error');
                    errorOccurred = true;
                    return;
                    
                end
                %}
            end
        end
        
        %crossTime{s}.mcp
        %crossTime{s}.gamma
        %{
        % Plot bad pulses
        noGood(:,i)
        timeOfFlight(s)
        %assignin('base','pulse',pulse);
        
        if noGood(c,i) == 1
            'entering if'
            c
            if c == 1
                assignin('base','pulseIn',pulse{c}(s).y);
                %assignin('base','foundFraction',foundFraction);
                %evalin('base','figure(ax)');
                evalin('base','plot(ax2,pulseIn);');
                
            elseif c == 2
                
                assignin('base','pulseIn',pulse{c}(s).y);
                %assignin('base','foundFraction',foundFraction);
                %evalin('base','figure(ax)');
                evalin('base','plot(ax3,pulseIn);');
            end
            %elseif (timeOfFlight(s) == 1 || timeOfFlight(s) == 0 || noGood(2,i) == 1) && c == 2
            %   assignin('base','pulseIn',pulse{c}(s).y);
            %assignin('base','foundFraction',foundFraction);
            %evalin('base','figure(ax)');
            %evalin('base','plot(ax3,pulseIn);');
            
        end
        %}
    end
    
    % Finding time difference between fractional value of MCP pulse and
    % x-crossover value of the gamma pulse.
    
    %try
        
        if sum(TypeIn(:) == 1) ~= 0
            
            if multiStopIn == '1'
                numPeaks = 2
                timeOfFlight{s} = NaN(1,numPeaks);
                'NO GOOD'
                noGood(:,s)
                noGoodS(:,s)
                if sum(noGood(:,s)) == 0 && sum(noGoodS(:,s)) == 0
                    'NUMPEAKS'
                    
                    for p=1:numPeaks
                        'ToF CALCULATION'
                        timeOfFlight{s}(p) = (crossTime{s}.mcp(p) - crossTime{s}.gamma(1))
                        
                        if ANN == 0 || (ANN == 1 && isempty(activeNeuron) == 0 && activeNeuron ~= 0)
                            'SUBTRACTING MULTI';
                            timeOfFlight{s}(p) = timeOfFlight{s}(p);% - ELETMatrix(activeNeuron,3);
                            
                        else
                            timeOfFlight(s) = NaN;
                            'multi-stop neuron error';
                            
                        end
                        
                        if timeOfFlight{s}(p) < -250*10.^(-9)
                            timeOfFlight{s}(p) = NaN;
                            'multi-stop out of range';
                            
                        end
                        
                    end
                    
                else
                    timeOfFlight{s}(1) = NaN;
                    'multi-stop no good';
                    
                end
                
            else
                noGood(:,s)
                noGoodS(:,s)
                if sum(noGood(:,s)) == 0 && sum(noGoodS(:,s)) == 0
                    'prior calculation';
                    timeOfFlight(s) = (crossTime{s}.mcp(1) - crossTime{s}.gamma(1));
                    %{
                    'mcp'
                    crossTime{s}.mcp(1)
                    'gamma'
                    crossTime{s}.gamma(1)
                    %}
                    if ANN == 0 || (ANN == 1 && isempty(activeNeuron) == 0 && activeNeuron ~= 0)
                        'SUBTRACTING';
                        timeOfFlight(s) = timeOfFlight(s);% - best(3,activeNeuron);%ELETMatrix(activeNeuron,3);
                        
                        
                        try
                                assignin('base','activeNeuron',activeNeuron);
                                assignin('base','timeOfFlight',timeOfFlight(s));
                                evalin('base','index = find(allNeurons.neurons{activeNeuron}(:,numberRuns) == 0);');
                                evalin('base','allNeurons.neurons{activeNeuron}(index(1),numberRuns) = timeOfFlight;');
                                
                        catch
                            'vector full'
                            
                        end
                        
                    else
                        timeOfFlight(s) = NaN;
                        'NOT SUBTRACTING'
                    end
                    
                else
                    timeOfFlight(s) = NaN;
                    'noGood'
                    
                end
                
                if timeOfFlight(s) < -250*10.^(-9)
                    timeOfFlight(s)
                    timeOfFlight(s) = NaN;
                    'out of range'
                    
                end
                
            end
            
        elseif sum(TypeIn(:) == 3) == 1 && sum(timingTypeIn(:) ~= 3) == 1 && 1==0
            
            if sum(noGood(:,s)) == 0 && sum(noGoodS(:,s)) == 0
                %timeOfFlight(s) = (crossTime{s}.gamma1(1) - crossTime{s}.gamma2(1));
                %figure
                inputIndex = 18750;
                %[~,inputIndex] = min(abs(pulse{1}(s).x-crossTime{s}.gamma1(1)));
                [output,value,maximum] = neuralNetworkFiltering(pulse{1}(s),inputIndex);
                
                if abs(crossTime{s}.gamma1(1)) < 400E-9 && ...
                        maximum ~= 0
                    
                        assignin('base','output',output);
                        assignin('base','value',value);
                        assignin('base','maximum',maximum);
                        evalin('base','inputVectors = horzcat(inputVectors,output);');
                        evalin('base','targetValues = horzcat(targetValues,[value;maximum]);');
                        
                        plot(pulse{1}(s).y,'-p','MarkerIndices',inputIndex, ...
                            'MarkerFaceColor','red','MarkerSize',15)
                end
                
            else
                timeOfFlight(s) = NaN;
                
            end
            
            if timeOfFlight(s) < -250*10.^(-9)
                timeOfFlight(s) = NaN;
                
            end
            
        else
            timeOfFlight(s) = NaN;
            
        end
       %{ 
    catch
        errordlg('Check channel timing type.','Timing Error');
        errorOccurred = true;
        return;
        
    end
    %}
    %clearvars pulse;
    s = s + 1;
    pause(0.0001);
    
end

%baselineAverage = baselineAverage/goodCounter;

end