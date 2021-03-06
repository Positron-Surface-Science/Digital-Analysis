function [VoutMax,timeOfFlight,errorOccurred,noFilesFound,dataFilter,baselineAverage] = ...
    tripleCoincidence_GeTiming_Trapezoidal_Functional_Input(selectPathIn,...
    iIn,analysisTypeIn,multiStopIn,multiStopConditionsIn,TypeIn,timingTypeIn,shapingTypeIn,...
    FFTIn,coinCh,ggc,RTFLIn,RTFUpIn,TFTopIn,TFRFIn,TFTopInFast,TFRFInFast,BPFSwitch,dataFilter,...
    DiscriminationVector,ParameterMatrix,ClusteringNetwork,neuron,ELETParam,ELETUpper,ELETLower)

%---------------------------------------------------------------------------
% Primary coincidence function. Reads the waveforms from the Lecroy
% oscilloscope files, then sends them to paes and trapezoidalfilter to
% calculate the time of flight and gamma values.
%---------------------------------------------------------------------------

selectPath = selectPathIn;%'J:\Lecroy\Coi';%uigetdir;

numChannels = 4;

%digits(2)
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
if isstruct(ParameterMatrix)
    a = fieldnames(ParameterMatrix);
    ParameterMatrix = getfield(ParameterMatrix,a{1});
end
if isempty(ggc)
    ggc = 0;
end
%ggc = 1;
exact2 = 1;
timingTypeChanged = false;
numPeaks = 0;
t = [];
k = [];

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
            if (timingTypeIn(c) == 4)
                ANN(c) = 1;
                %ggc = 1;
                %timingTypeIn(c) = 0;
                timingTypeChanged = true;
                'ANN ACTIVE FOR TIMING';
                
            elseif TypeIn(c) == 3 && shapingTypeIn(c) == 0 && timingTypeIn(c) ~= 0
                ANN(c) = 0;
                'ANN ACTIVE FOR SHAPING';
                
            else
                ANN(c) = 0;
                %activeNeuron = 0;
                'ANN NOT ACTIVE';
                ELETMatrix = [0 0];
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
                %try
                notThisOne = false;
                    try
                        lo = dir([selectPath,'\*00001.trc']);
                        k = strfind(lo(1).name,'Trace');
                        zo = [];
                        t = [];
                    catch
                       notThisOne = true; 
                       
                    end
                    
                    if notThisOne == true || isempty(k)%(true)
                        
                        try
                            lo = [];
                            k = [];
                            zo = dir([selectPath,'\*pulse_1.txt']);
                            t = strfind(zo(1).name,'pulse_');
                            noT = false;
                            
                        catch
                        end
                        
                        if isempty(t)
                            
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
                            
                        elseif ~isempty(t)
                            file = ['cebr3_vs_plastic_pulse_pulse_',stringconversion(n)];
                            
                        end
                        
                    end
                    
                    if ~isempty(k)%(false)
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
                    
                    if isempty(t)
                        fileName = [selectPath,'\',channelNumber{c},appendage,'.trc'];
                        
                    elseif ~isempty(t)
                        fileName = [selectPath,'\',file,'.txt'];
                        
                    end
                    
                %catch
                   % noFile = true;
                %end
                
                % try to read the file; if it doesn't exist, try again twenty
                % times every five seconds, then break and try to read the
                % next file if its existence never occurs.
                
                tCounter = 1;
                
                while noFile == true && noFileCounter <= 7 && tCounter <= 10
                    'No File Loop';
                    try
                        if isempty(t)
                            pulse{c}(s) = waveform.ReadLeCroyBinaryWaveform(fileName);
                            %pulse{c}(s).y = single(pulse{c}(s).y);
                            %pulse{c}(s).y
                            noFile = false;
                            
                        elseif ~isempty(t)
                            
                            allPulses = importdata(fileName, '\t', 1);
                            assignin('base','allPulses',allPulses);
                            
                            if c == 1
                                pulse{c}(s).x = allPulses.data(:,1)*10.^(-9);
                                pulse{c}(s).y = allPulses.data(:,2)*10.^(-3);
                                noFile = false;
                                
                            elseif c == 2
                                pulse{c}(s).x = allPulses.data(:,3)*10.^(-9);
                                pulse{c}(s).y = allPulses.data(:,4)*10.^(-3);
                                noFile = false;
                                
                            end
                        end
                        
                    catch
                        %errordlg('No .trc files within folder.','Incorrect Folder.');
                        noFile = true;
                        noFileCounter = noFileCounter + 1;
                        pause(1);
                        
                    end
                    
                    if noFileCounter >= 5 || tCounter >= 10
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
                    
                    tCounter = tCounter + 1;
                    
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
                
                if isempty(t)
                    offset = -pulse{c}(s).info.OFFSET;
                    Fs = pulse{c}(s).desc.fs;
                    T = (1/Fs);
                    oldPulse = pulse{c}(s).y;
                    
                elseif ~isempty(t)
                    offset = 13E-9;
                    Fs = 5E9;
                    T = (1/Fs)
                    oldPulse = pulse{c}(s).y;
                    
                end
                
            end
            
            if (false)
                [pulse{c}(s).x pulse{c}(s).y,tCell] = DRS4BinaryDataOpenTest(c,n);
                
                offset = tCell;
                Fs = 1.224;
                T = (1/Fs);
                pulse{c}(s).desc.Ts = T;
                oldPulse = pulse{c}(s).y;
                
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
                    try
                        baseline = mean(pulse{c}(s).y(round((offset-999E-9)/T):round((offset-500E-9)/T)));
                        pulse{c}(s).y = pulse{c}(s).y - baseline;
                        'BASELINE CORRECTION PERFORMED'
                        
                    catch
                        
                    end
                    
                end
                
            catch
                'bpf error';
            end
            
            pulse1 = single(pulse{c}(s).y(7301:8200));
            
            %if max(pulse1) <= 0.0595
                %assignin('base','pulse1',pulse1);
                %evalin('base','pulsesBack=horzcat(pulsesBack,pulse1);');
                
            %end
            
            if c == 2 && false
                %plot(pulse{c}(s).y(45:500))
                %assignin('base','trace',pulse{c}(s).y(45:500));
                %evalin('base','traces = horzcat(traces,trace);');
                %assignin('base','p',pulse{c}(s).y(18626:19063)/max(pulse{c}(s).y(18626:19063)));
                %assignin('base','xTrace',pulse{c}(s).x(18626:19063));
                netTrace = pulse{c}(s).y(round(10*10.^(-6)/T):round(35*10.^(-6)/T));
                assignin('base','netTrace',netTrace);
                %evalin('base','netTraces = horzcat(netTraces,netTrace);');
            end
            %{
            if c == 2
                netTrace2 = pulse{2}(s).y(45:500);
                
            end
            %}
            %{
            if TypeIn(1) == 4 && TypeIn(c) == 3
                assignin('base','trace',pulse{c}(s).y(round((offset-999E-9)/T):round((offset+500E-9)/T)));
                evalin('base','NaITraces = horzcat(NaITraces,trace(2550:3000));');
            
            end
            %}
            %{
            if TypeIn(c) == 1 && ggc == 1
                try
                [m,i] = max(pulse{c}(s).y);
                assignin('base','trace',pulse{c}(s).y(i-100:i+100));
                mcpMatrix = evalin('base','net(trace);');
                exact1 = evalin('base','exact');
                'EXACT 2'
                exact2 = exact1*mcpMatrix
                catch
                    exact2 = 0;
                end
            end
            %}
            
            
            %---------------------------------------------------------------------------------------
            % ANN Analysis.
            %---------------------------------------------------------------------------------------
            
            %try
            if ANN(c) == 1
                'ANN ON'
                p = 1;
                q = floor((1.25*p)/0.5);
                FsN = (p/q)*Fs;
                TsN = 1/FsN;
                %{
                try
                    baseline = mean(oldPulse(round((offset-999E-9)/T):round((offset-500E-9)/T)));
                    oldPulse = oldPulse - baseline;
                catch
                    'bpf error 2'
                end
                %}
                %if evalin('base','flags;') == 1
                
                %assignin('base','trace',pulse{c}(s).y(2500:4000));
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
                try
                    %trace2 = pulse{c}(s).y(1:5000);
                    %assignin('base','trace2',trace2);
                    netTrace = pulse{c}(s).y(45:500);
                    %assignin('base','netTrace',netTrace);
                    %netTrace = netTrace(10:55);
                    %plot(netTrace)
                    %netTrace = pulse{c}(s).y(18626:19063); %pulse{c}(s).y(2500:4000);
                    %netTracex = pulse{c}(s).x(2500:4000);
                    %assignin(
                    %oldPulse(18626:19063);%
                    %netTrace = netTrace/max(netTrace);
                    %'HEREHERHERHEHER'
                    %shapingTypeIn(c)
                    if shapingTypeIn(c) == 0
                        'HPGe SHAPING'
                        netTrace = pulse{c}(s).y(round(10*10.^(-6)/T):round(35*10.^(-6)/T));
                    
                        netTrace = netTrace(501:1000);
                        
                        assignin('base','netTrace',netTrace);
                        %evalin('base','pulses = horzcat(pulses,netTrace)');
                        
                    end
                    
                catch
                    netTrace = [];
                    'no trace'
                    
                end
                
                %assignin('base','xTrace',pulse{c}(s).x(18626:19063));
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
                %network(net);
                %best = ParameterMatrix;
                %net = ClusteringNetwork; %evalin('base','net;'); %
                net = evalin('base','net');
                %a = fieldnames(DiscriminationVector);
                %exact = getfield(DiscriminationVector,a{1});
                %mean(single(netTrace*1E4))
                %assignin('base', 'paramTest', best);
                %plot(netTrace)
                out = net(single(netTrace)); %net(single(netTrace*1E4));
                %'elet2'
                %net2 = evalin('base','net4;');
                %p2 = pulse{c}(s).y(18626:19063)/max(pulse{c}(s).y(18626:19063));
                %netTrace2 = p2*1E2;
                %elet2 = predict(net2,netTrace2);
                %'ELET 2'
                %elet2(3) = -elet2(3)*1E-7;
                %elet2
                'Active Neuron'
                activeNeuron = find(out == 1)
                
                
                %nnz
                
                
                %pause(1)
                %assignin('base','netTrace',netTrace);
                %evalin('base','plot([netTrace,tracesSingle(:,1)]);');
                %'timing type'
                %timingTypeIn(c)
                
                if timingTypeIn(c) ~= 3
                    %net = evalin('base','net;');
                    %{
                    if ~ismember(activeNeuron, ...
                            [633 634 635 636 637 638 ...
                            664 665 666 667 668 669 670 671 672 ...
                            698 702 703 704 ...
                            733 734 735 736 ...
                            765 766 767 768 ...
                            799 800 ...
                            830 831 832 ...
                            1022 1023 1024 ...
                            990 991 992 ...
                            958 959 960 ...
                            926 927 928 ...
                            894 895 896 ...
                            862 863 864])
                    %}
                            %{
                            [1007 1008 1009 1010 1011 1012 1013 1014 1015 ...
                            974 975 976 977 979 980 981 982 983 ...
                            942 943 944 946 947 948 949 950 951 ...
                            909 910 911 912 913 914 915 916 917 918 919 ...
                            876 877 878 879 880 881 882 883 884 885 886 ...
                            844 845 846 847 848 849 850 851 852 ...
                            812 813 814 815 816 817 818 819 ...
                            783 784 785 786])
                            %}
                            %{
                            Upper Right
                            [633 634 635 636 637 638 ...
                            664 665 666 667 668 669 670 671 672 ...
                            698 702 703 704 ...
                            733 734 735 736 ...
                            765 766 767 768 ...
                            799 800 ...
                            830 831 832 ...
                            1022 1023 1024 ...
                            990 991 992 ...
                            958 959 960 ...
                            926 927 928 ...
                            894 895 896 ...
                            862 863 864])
                            %}
                        
                            %%% Lower Left
                            %[3 4 5 6 7 8 ...
                            %34 35 36 37 38 39 ...
                            %65 66 67 68 69 70 71 ...
                            %97 98 99 100 101])
                            
                            %activeNeuron = 0;
                            
                    %end
                    %}
                elseif timingTypeIn(c) == 3
                    %net = evalin('base','net;');%getfield(ClusteringNetwork,a{1});
                    %{
                    'IS MEMBER?'
                    
                    if ~ismember(activeNeuron, ...
                            [1007 1008 1009 1010 1011 1012 1013 1014 1015 ...
                            974 975 976 977 979 980 981 982 983 ...
                            942 943 944 946 947 948 949 950 951 ...
                            909 910 911 912 913 914 915 916 917 918 919 ...
                            876 877 878 879 880 881 882 883 884 885 886 ...
                            844 845 846 847 848 849 850 851 852 ...
                            812 813 814 815 816 817 818 819 ...
                            783 784 785 786])
                            %{
                            Upper Right
                            [633 634 635 636 637 638 ...
                            664 665 666 667 668 669 670 671 672 ...
                            698 702 703 704 ...
                            733 734 735 736 ...
                            765 766 767 768 ...
                            799 800 ...
                            830 831 832 ...
                            1022 1023 1024 ...
                            990 991 992 ...
                            958 959 960 ...
                            926 927 928 ...
                            894 895 896 ...
                            862 863 864])
                            %}
                        
                            %%% Lower Left
                            %[3 4 5 6 7 8 ...
                            %34 35 36 37 38 39 ...
                            %65 66 67 68 69 70 71 ...
                            %97 98 99 100 101])
                            
                        activeNeuron = 0;
                        
                    end
                    %}
                end
                
                %timing = evalin('base','net3(activeNeuron)');
                assignin('base','activeNeuron',activeNeuron);
                
                %if timingTypeIn(c) ~= 3 && activeNeuron ~= 0
                
                %best = evalin('base','best;');
                
                %'ELET PARAMETERS'
                %ELETParam = best(1:2,activeNeuron);
                %ParameterMatrix = evalin('base','best2;');
                %ELETParam
                'ELET PARAM';
                ELETMatrix = ParameterMatrix(1:2,activeNeuron);
                ELETMatrix
                %best(1:2,activeNeuron);
                
                %if best(1,activeNeuron) == 0
                %ELETMatrix(activeNeuron,1:2) = [7 9];
                
                %end
                
                'FWHM';
                %ELETMatrix(activeNeuron,4)
                
                %elseif timingTypeIn(c) ~= 3
                
                %activeNeuron = 1;
                %ELETMatrix(activeNeuron,1:2) = [0 0];
                
                %end
                
                %'exact'
                %exact(activeNeuron)
                
                %ELETMatrix(activeNeuron,1:2) = best(activeNeuron,1:2);
                %{
                if  ggc == 1 && (isempty(activeNeuron) || ...
                        isempty(neuron) || activeNeuron ~= neuron)
                    %activeNeuron = 0;
                    
                elseif ggc == 1 && isempty(activeNeuron) == 0 && ...
                        activeNeuron == neuron
                    %ELETMatrix(activeNeuron,1:2) = ELETParam;
                    %ELETMatrix(activeNeuron,3) = 0;
                    'EQUAL'
                    
                end
                %}
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
                %activeNeuron = 0;
                
            %end
            
            %---------------------------------------------------------------------------------------
            % HPGe pulse analysis (trapezoidal).
            %---------------------------------------------------------------------------------------
            
            if (TypeIn(c) == 3 || TypeIn(c) == 2) && shapingTypeIn(c) ~= 8 && ...
                    ((sum(ANN) == 1 && activeNeuron ~= 0) || sum(ANN) == 0)
                
                [TPA,~,noGoodS(c,s)] = trapezoidalFilter(pulse{c}(s),oldPulse,shapingTypeIn(c),TFRFIn(c),TFTopIn(c),0);
                
                VoutMax{c}(s) = TPA;
                
                if noGoodS(c,s) == false
                    goodCounter = goodCounter + 1;
                    %baselineAverage = baselineAverage + baseline;
                    
                end
                
            else
                VoutMax{c}(s) = 0;
                noGood(c,s) = 1;
                
            end
            
            %---------------------------------------------------------------------------------------
            % PAES analysis.
            %---------------------------------------------------------------------------------------
            
            if timingTypeIn(c) ~= 3 && (c == coinCh(1) || c == coinCh(2))
                
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
                %plot(pulse{c}(s).y)
                
                if TypeIn(c) == 3 && timingTypeIn(c) == 2
                    [~,vOut] = trapezoidalFilter(pulse{c}(s),oldPulse,1,2E-6,0.1E-6,1);
                    %plot(pulse{c}(s).y);
                    try
                        plot(vOut(15000:20000,:));
                        pulse{c}(s).y = vOut(:,2);
                    catch
                        pulse{c}(s).y = [];
                        'no pulse'
                    end
                    
                end
                
                %plot(pulse{c}(s).x,pulse{c}(s).y)
                %pause(1)
                [crossTimeOut,noGood(c,s),numPeaks,~] = paes_unedited(timingTypeIn(c),TypeIn(c),c,multiStopIn,...
                    multiStopConditionsIn,numPeaks,pulse{c}(s),RTFLIn,RTFUpIn,ELETMatrix, ...
                    DiscriminationVector,neuron,activeNeuron,ANN,ggc,ELETUpper,ELETLower,t);
                
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
        %{
        if TypeIn(c) == 1
            
            assignin('base','pulse_BaF2',pulse{c}(s).y(18626:19063)/max(pulse{c}(s).y(18626:19063)));
            assignin('base','time_BaF2',crossTimeOut);
            evalin('base','pulses_BaF2 = horzcat(pulses_BaF2, pulse_BaF2);');
            evalin('base','times_BaF2 = horzcat(times_BaF2, time_BaF2);');
            
        elseif TypeIn(c) == 3
            
            assignin('base','pulse_HPGe',pulse{c}(s).y(18626:19063)/max(pulse{c}(s).y(18626:19063)));
            assignin('base','neuron_HPGe',activeNeuron);
            evalin('base','pulses_HPGe = horzcat(pulses_HPGe, pulse_HPGe);');
            evalin('base','neurons_HPGe = horzcat(neurons_HPGe, neuron_HPGe);');
            
        end
        %}
    end
    
    % Finding time difference between fractional value of MCP pulse and
    % x-crossover value of the gamma pulse.
    
    %try
    
        if sum(TypeIn(:) == 1) ~= 0
            
            if multiStopIn == '1'
                %numPeaks = 2
                %{
                if multiStopConditionsIn == 6 || multiStopConditionsIn == 2
                    numPeaks = 2;
                    
                elseif multiStopConditionsIn == 7 || multiStopConditionsIn == 3
                    numPeaks = 3;
                    
                elseif multiStopConditionsIn == 1
                    numPeaks = 1;
                    
                end
                %}
                numPeaks = numel(crossTime{s}.mcp)
                
                timeOfFlight{s} = NaN(1,numPeaks);
                'NO GOOD';
                noGood(:,s);
                noGoodS(:,s);
                if sum(noGood(:,s)) == 0 && sum(noGoodS(:,s)) == 0
                    
                    for p=1:numPeaks
                        'ToF CALCULATION';
                        timeOfFlight{s}(p) = (crossTime{s}.mcp(p) - crossTime{s}.gamma(1));
                        
                        if numPeaks == 3
                            %'NO GOOD'
                            %noGood(:,s)
                            %noGoodS(:,s)
                            timeOfFlight{s}(p)
                            %'mcp 1'
                            %crossTime{s}.mcp(p)
                            %'gamma 1'
                            %crossTime{s}.gamma(1)
                            
                            %pause(2)
                            
                        end
                        
                        'mcp 1';
                        crossTime{s}.mcp(p);
                        'gamma 1';
                        crossTime{s}.gamma(1);
                        
                        if (sum(ANN) == 1 && isempty(activeNeuron) == 0 && activeNeuron ~= 0)
                            'SUBTRACTING MULTI'
                            timeOfFlight{s}(p) = timeOfFlight{s}(p) - ParameterMatrix(3,activeNeuron);
                            
                        elseif sum(ANN) ~= 0 && (isempty(activeNeuron) || activeNeuron == 0)
                            timeOfFlight(s) = NaN;
                            'multi-stop neuron error';
                            
                        end
                        
                        if timeOfFlight{s}(p) < -250*10.^(-9)
                            timeOfFlight{s}(p) = NaN;
                            'multi-stop out of range'
                            
                        end
                        
                    end
                    
                    try
                        if ~isnan(timeOfFlight{s})
                            assignin('base','timeOfFlight',timeOfFlight{s});
                            evalin('base','ToF = vertcat(ToF,timeOfFlight);');
                            
                        end
                    catch
                    end
                    
                else
                    timeOfFlight{s}(1) = NaN;
                    'multi-stop no good';
                    
                end
                
            else
                try
                noGood(:,s)
                noGoodS(:,s)
                if sum(noGood(:,s)) == 0 && sum(noGoodS(:,s)) == 0 && exact2 == 1
                    'prior calculation'
                    timeOfFlight(s) = (crossTime{s}.mcp(1) - crossTime{s}.gamma(1))
                    
                    'mcp 2';
                    crossTime{s}.mcp(1);
                    'gamma 2';
                    crossTime{s}.gamma(1);
                    %crossTime{s}.gamma(1)
                    %
                    %sum(ANN)
                    %isempty(activeNeuron)
                    %activeNeuron
                    
                    if (sum(ANN) == 1 && ~isempty(activeNeuron) && activeNeuron ~= 0)% && ...
                            %VoutMax{2}(s) >= 5.94 && VoutMax{2}(s) <= 7.23%6.713%best(4,activeNeuron) <= 10E-9 && VoutMax{2}(s) >= 6.2 && VoutMax{2}(s) <= 7.0 6.2
                            
                        'SUBTRACTING'
                        timeOfFlight(s) = timeOfFlight(s) - ParameterMatrix(3,activeNeuron) %ELETMatrix(activeNeuron,3);
                        %{
                        assignin('base','netTrace',netTrace);
                        assignin('base','netTrace2',netTrace2);
                        
                        if timeOfFlight(s) <= -1E-10
                            evalin('base','lowPeak = horzcat(lowPeak,netTrace);');
                            evalin('base','lowPeak2 = horzcat(lowPeak2,netTrace2);');
                            
                        elseif timeOfFlight(s) > -1E-10 && timeOfFlight(s) < 1E-10
                            evalin('base','midPeak = horzcat(midPeak,netTrace);');
                            evalin('base','midPeak2 = horzcat(midPeak2,netTrace2);');
                            
                        elseif timeOfFlight(s) >= 1E-10
                            evalin('base','highPeak = horzcat(highPeak,netTrace);');
                            evalin('base','highPeak2 = horzcat(highPeak2,netTrace2);');
                            
                        end
                        %}
                        if true
                            %try
                            'Assigning ToF to activeNeuron'
                            activeNeuron
                            timeOfFlight(s)
                                assignin('base','activeNeuron',activeNeuron);
                                assignin('base','timeOfFlight',timeOfFlight(s));
                                evalin('base','index = find(allNeurons.neurons{activeNeuron}(:,numberRuns) == 0);');
                                evalin('base','allNeurons.neurons{activeNeuron}(index(1),numberRuns) = timeOfFlight;');
                                
                            %catch
                            %    'vector full'
                                
                            %end
                        end
                    
                    elseif ANN ~= 0 && (isempty(activeNeuron) || activeNeuron == 0)
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
                catch
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
    %{
    if abs(timeOfFlight(s)) <= 3E-9
        evalin('base','pulses = horzcat(pulses,p);');
        evalin('base','elets = horzcat(elets, best3(1:3,activeNeuron));');
        
    end
    %}
    s = s + 1;
    pause(0.0001);
    
end

%baselineAverage = baselineAverage/goodCounter;

end