function [crossTime,noGood,numberOfPeaks,riseTime] = paes_unedited(timingIn,TypeIn,c,...
    multiStopIn,multiStopConditionsIn,numPeaks,pulseIn,RTFLIn,RTFUpIn,ELETMatrix,...
    DiscriminationVector,neuron,activeNeuron,ANN,ggc,ELETUpper,ELETLower,t)

% --------------------------------------------------------------------------
% 'PAES' function which calculates the crossing "pick-off" time of detector
% pulses.
% --------------------------------------------------------------------------
%pulseIn.y = single(pulseIn.y);
%pulseIn.x = single(pulseIn.x);
pulseIny = pulseIn.y;
pulseInx = pulseIn.x;
tempOn = false;

if ~isempty(t) && TypeIn == 1
    TypeIn = 3;
    tempOn = true;
    
end

try
    dN = pulseIn.desc.Ts; % Time per channel
    
catch
    dN = 1/(5E9);
    
end

riseTime = 1;
foundFraction = [];
foundFraction1 = [];
foundFraction2 = [];
count = 0;
count1 = 0;
count2 = 0;
foundFrac = 1;
foundFrac1 = 1;
foundFrac2 = 1;
goodRiseTime = false;
findingNoise = false;
neuronTemp = neuron;
scintPulse = 0;
gggc = 0;
numberOfPeaks = 1;
rejecting = false;

%ELETMatrix = ELETMatrix';
%a = fieldnames(DiscriminationVector);
%exact = getfield(DiscriminationVector,a{1});
% ----------------------------------
%exact(exact == 0) = 1;
% ----------------------------------

if TypeIn == 1 && ggc == 0
    neuronTemp = 0;
    
    % Considering if a pulse is negative- or positive-going.
    if max(pulseIny) < abs(min(pulseIny))
        pulseIny = -pulseIny;
        
    end
    
    lastwarn('');
    % Set to detect peaks 15 ns apart (15*10.^(-9)/dN)
    [VMin{c},VMinIndex{c}] = findpeaks(pulseIny(50:numel(pulseIny)-50),...
        'MinPeakDistance',round((14*10.^(-9))/dN), ...
        'MinPeakHeight',1*10.^(-3),'MinPeakProminence',0.25*10.^(-3), ...
        'MinPeakWidth',round((10*10.^(-9))/dN),'WidthReference', ...
        'halfheight','MaxPeakWidth',round((75*10.^(-9)/dN)));
    
    [warnMsg, ~] = lastwarn;
    VMinIndex{c} = VMinIndex{c} + 50;
    numPeaks = numel(VMinIndex{c});
    numberOfPeaks = numPeaks;
    crossTime = NaN(numPeaks,1);
    'numPeaks';
    numPeaks;
    
    %plot(pulseIny(VMinIndex{c}-50:VMinIndex{c}+50))
    
    if numPeaks == 0 && multiStopConditionsIn ~= 6  && ...
            multiStopConditionsIn ~= 7
        
        [VMin{c},VMinIndex{c}] = max(pulseIny(50:numel(pulseIny)-50));
        VMinIndex{c} = VMinIndex{c} + 50;
        warnMsg = [];
        numPeaks = 1;
        numberOfPeaks = numPeaks;
        crossTime = NaN(numPeaks,1);
        'weird error'
        
    end
    
    if multiStopIn == '1' && multiStopConditionsIn == 7 && numPeaks == 2
        %multiStopConditionsIn = 3;
        
    elseif multiStopIn == '1' && multiStopConditionsIn == 7 && numPeaks == 3
        multiStopConditionsIn = 3;
        
    elseif multiStopIn == '1' && multiStopConditionsIn == 6
        multiStopConditionsIn = 2;
        
    end
   
    if multiStopIn == '0'
        numPeaks = 1;
        numberOfPeaks = 1;
        
    end
    
    %{
     if numPeaks == 3
        TypeIn == 1
        multiStopIn == '1'
        numel(warnMsg) > 0
        multiStopConditionsIn ~= 4
        multiStopConditionsIn ~= 5
        'NUMPEAKS'
        multiStopConditionsIn
        numPeaks
        multiStopConditionsIn ~= numPeaks
        (multiStopConditionsIn == 5 && (numPeaks == 1 || numPeaks == 0))
        numPeaks > 5
        
        rejecting = true;
        
        pause(15)
        
    end
    %}
    %if numPeaks == 0
        %figure
     %   plot(pulseIn.y)
    %end
    
elseif TypeIn == 1 && ggc == 1
    neuronTemp = 0;
    [VMin{c},VMinIndex{c}] = max(pulseIny(50:numel(pulseIny)-50));
    [VMin2{c},VMinIndex2{c}] = min(pulseIny(50:numel(pulseIny)-50));
    
    % Considering if a pulse is negative- or positive-going.
    if abs(VMin2{c}) > abs(VMin{c})
        VMin{c} = abs(VMin2{c});
        VMinIndex{c} = VMinIndex2{c};
        pulseIny = abs(pulseIny);
        
    end
    
    VMinIndex{c} = VMinIndex{c} + 50;
    warnMsg = [];
    crossTime = NaN;
    numPeaks = 1;
    %plot(pulseIny)
    
elseif TypeIn ~= 1
    'VMin';
       %{         
    [VMin{c},VMinIndex{c}] = max(earlyOne(10:numel(earlyOne)-50));
    [VMin2{c},VMinIndex2{c}] = min(earlyOne(10:numel()-50));
         %}
               
    [VMin{c},VMinIndex{c}] = max(pulseIny(10:numel(pulseIny)-50));
    [VMin2{c},VMinIndex2{c}] = min(pulseIny(10:numel(pulseIny)-50));
    
    %pause(2)
    
    % Considering if a pulse is negative- or positive-going.
    if abs(VMin2{c}) > abs(VMin{c})
        VMin{c} = abs(VMin2{c});
        VMinIndex{c} = VMinIndex2{c};
        pulseIny = abs(pulseIny);
        scintPulse = 1;
        
    end
    
    %VMin{c} = (VMin{c} + mean(pulseIny(VMinIndex{c}-3:VMinIndex{c}+3)))/2;
    VMinIndex{c} = VMinIndex{c} + 10;
    warnMsg = [];
    crossTime = NaN;
    numPeaks = 1;
    
    %{
elseif TypeIn == 1
    [VMin{c},VMinIndex{c}] = max(pulseIny(50:numel(pulseIny)-50));
    
    VMinIndex{c} = VMinIndex{c} + 50;
    warnMsg = [];
    crossTime = 1;
    %}
end

noGood = false;
VMinFraction = zeros(1,numPeaks);
VMinFraction1 = zeros(1,numPeaks);
VMinFraction2 = zeros(1,numPeaks);

numberOfPeaks = numPeaks;

for s=1:numberOfPeaks
        %isempty(VMinIndex{c})
    if isempty(activeNeuron) ||  ...
            isempty(VMinIndex{c}) || ...
            VMinIndex{c}(1) - 10 <= 0
        crossTime(s) = NaN;
        numPeaks = 0;
        'neuron error 1'
        break;
    
    end
        
    if TypeIn == 3 && ANN(c) == 1 && activeNeuron == 0
        crossTime(s) = NaN;
        numPeaks = 0;
        'neuron error 2'
        break;
        
    end
    
    % MCP ANN
    if TypeIn == 1
        try
            %{
            assignin('base','pulseIny',pulseIny(VMinIndex{c}-100:VMinIndex{c}+100));
            out = evalin('base','net(pulseIny);');
            in = evalin('base','next;');
            if in*out == 0
                crossTime(s) = NaN;
                numPeaks = numPeaks - 1;
                riseTime = 1;
                break;
            end
            %}
        catch
            '?'
        end
    end
    
    %'s'
    %s
    % Integrating peak and region before peak to reduce periodic noise
    % pulses.
    if (TypeIn == 1 && ggc == 0) && VMinIndex{c}(s)-round((50*10.^(-9)/dN)) > 0
        
        integration = trapz(pulseIny(VMinIndex{c}(s)-round((50*10.^(-9)/dN)) ...
            :VMinIndex{c}(s)));
        
    elseif TypeIn == 1
        integration = 0;
        %'TypeIn 1'
        
    end
    
    % WAS integration < 0.005
    if (TypeIn == 1 && ggc == 0) && (integration < 0.005) % || integration >= 0.01
       numPeaks = numPeaks - 1;
       crossTime(s) = NaN;
       'integration'
       integration
       continue;
        
    end
    
    % Preventing multiple analyses on gamma pulses
    if TypeIn ~= 1 && s > 1
        's > 1'
        break;
        
    else
        % Fraction for CFD, ELET, and RTF.
        VMinFraction = 0.01*(ELETLower(c)*VMin{c}(s));
        VMinFraction1 = 0.01*(ELETLower(c)*VMin{c}(s));
        VMinFraction2 = 0.01*(ELETUpper(c)*VMin{c}(s));
        
        if ANN(c) == 1
            %assignin('base','ELETMatrix',ELETMatrix);
            VMinFraction1 = 0.01*(ELETMatrix(1)*VMin{c}(s));
            VMinFraction2 = 0.01*(ELETMatrix(2)*VMin{c}(s));
            
            if ELETMatrix(1) == ELETMatrix(2)
                cfdANN = true;
                VMinFraction = 0.01*(ELETMatrix(1)*VMin{c}(s));
                
            end
            
            %VMinFraction1 = (0.07*VMin{c}(s));
            %VMinFraction2 = (0.09*VMin{c}(s));
            
        end
        
        riseTimeLower = (0.1*VMin{c}(s));
        riseTimeHigher = (0.9*VMin{c}(s));
        %riseTimeLower = (0.1*VMin(c,s));
        %riseTimeHigher = (0.9*VMin(c,s));
        
        % Charging from the peak towards lower time to find the closest value to
        % the fractional voltage.
        
        %VLower2 = 0.001;
        %VUpper2 = 1000;
        V = VMin{c}(s);
        
        if TypeIn == 1
            VLower = 0.0007;
        else
            VLower = 0.0007;
        end
        % Ensuring pulseIn indices are greater than 1 and that the pulseIn is
        % within a particular voltage range.
        T = V > VLower;% && V < VUpper;
        %T2 = V > VLower2 && V < VUpper2;
        
    end
    
    % Preventing multiple analyses on gamma pulses.
    if TypeIn ~= 1 && s > 1
        'TypeIn & s'
        break;
        
    % ----------------------------------------
    % Constant fraction discrimination: finding channel corresponding to
    % fraction.
    % ----------------------------------------
    
    elseif  (timingIn == 1 || timingIn == 2 || (timingIn == 4 && cfdANN == true)) && T ~= 0 && VMinIndex{c}(s) > 1
        
        %{
        stop = 0;
        count = 0;
        
        % count cannot exceed 500 ns.
        while stop == 0 && VMinIndex{c}(s) > count && count < (500*10.^(-6)/pulseIn.desc.Ts)
            
            % Going down peak edge by 'count' until less than VMinFraction.
            if abs(pulseIny(VMinIndex{c}(s) - count)) <= abs(VMinFraction)
                foundFrac = (VMinIndex{c}(s) - count);
                stop = 1;
                
            else
                count = (count + 1);
                
            end
        end
        %}
        
        if TypeIn == 1 && ggc == 0 && scintPulse == 0 && tempOn == false
            'IF';
            if s >= 2 && pulseInx(VMinIndex{c}(s)) - pulseInx(VMinIndex{c}(s-1)) ...
                    <= 50*10.^(-9)
                
                [~, foundFrac] = min(abs(pulseIny(VMinIndex{c}(s-1) ...
                    :VMinIndex{c}(s)) - VMinFraction));
                
                foundFrac = foundFrac + VMinIndex{c}(s-1);
                foundFraction = horzcat(foundFraction,foundFrac);
                
            else
                
                [~, foundFrac] = min(abs(pulseIny(VMinIndex{c}(s)- ...
                    round(50*10.^(-9)/dN):VMinIndex{c}(s)) - VMinFraction));
                
                foundFrac = foundFrac + (VMinIndex{c}(s) - round(50*10.^(-9)/dN));
                foundFraction = horzcat(foundFraction,foundFrac);
                
            end
            
            %{
            [~, foundFrac] = min(abs(pulseIny(VMinIndex{c}(s-1)+ ...
                round(10*10.^(-9)/dN):VMinIndex{c}(s)) - VMinFraction));
                foundFrac = foundFrac + VMinIndex{c}(s-1) + round(20*10.^(-9)/dN);
            %}
            %assignin('base','pulseIn',pulseIny(VMinIndex{c}(s-1)+round(10*10.^(-9)/dN):numel(pulseIny)));
            
        else
            'ELSE';
            % WAS 700 ns
            
            %try
            
            earlyOne = smooth(pulseIny(1: ...
                VMinIndex{c}(s)), 5);
            
            leap = diff(smooth(diff(earlyOne), 5));
            
                [~, risingEdge] = max(leap(50:end-5));
                
                risingEdge = risingEdge + 45;
                
                %{
                [~, midPulse] = min(abs(pulseIny(risingEdge: ...
                    VMinIndex{c}(s)) - 0.25*VMin{c}(s)));
                
                midPulse = midPulse + risingEdge + 5;
                %}
                
                %size(leap)
                %size(pulseIny(1:VMinIndex{c}(s)))
                
                if tempOn == true
                    %plot([leap pulseIny(1:VMinIndex{c}(s)-2)])
                    
                end
                
                %{
                [~, foundFrac] = min(abs(pulseIny(VMinIndex{c}(s)-round(700E-9/dN): ...
                    VMinIndex{c}(s)) - VMinFraction));
                foundFrac = foundFrac + (VMinIndex{c}(s) - round(700E-9/dN));
                %foundFraction = horzcat(foundFraction,foundFrac);
                %}
                
                [~, foundFrac] = min(abs(pulseIny(risingEdge: ...
                    VMinIndex{c}(s)) - VMinFraction));
                
                if tempOn == true
                'FOUNDFRAC'
                c
                risingEdge
                foundFrac
                VMinIndex{c}(s)
                [~,maximum] = max(pulseIny);
                maximum
                end
                
                foundFrac = foundFrac + risingEdge;
                
             %{   
            catch
                
                [~, foundFrac] = min(abs(pulseIny(1: ...
                    VMinIndex{c}(s)) - VMinFraction));
                
            end
            %}
        end
        
    % ---------------------------------
    % ELET: finding channels corresponding to the two fractions.
    % ---------------------------------
    %activeNeuron
    elseif (timingIn == 0 || timingIn == 4) && T ~= 0 && ...
            VMinIndex{c}(s) > round(1E-9/dN)% && activeNeuron == 0% ...
            %|| exact(activeNeuron) == 1)
        
        %{
        stop1 = 0;
        stop2 = 0;
        
        while (stop1 == 0 || stop2 == 0) && VMinIndex{c}(s) > count2 && VMinIndex{c}(s) > count1 ...
                && count2 < 0.12*numel(pulseIny) && count1 < 0.12*numel(pulseIny)
            
            % Going down peak edge by 'count' until less than VMinFraction.
            if stop1 == 0 && abs(pulseIny(VMinIndex{c}(s) - count1)) <= abs(VMinFraction1)
                foundFrac1 = (VMinIndex{c}(s) - count1);
                stop1 = 1;
                
            elseif stop1 == 0
                count1 = (count1 + 1);
                
            end
            
            % Going down peak edge by 'count' until less than VMinFraction.
            if  stop2 == 0 && abs(pulseIny(VMinIndex{c}(s) - count2)) <= abs(VMinFraction2)
                foundFrac2 = (VMinIndex{c}(s) - count2);
                stop2 = 1;
                
            elseif stop2 == 0
                count2 = (count2 + 1);
                
            end
        end
        %}
        %pulse2 = pulseIny;
        %[pulseIny,error,~] = MLSmoothing(pulseIny,dN);
        %{
        if error > 1E-2
            noGood = true;
            riseTime = 1;
            numPeaks = 0;
            break;
            
        end
        %}
        
        try
            %{
            if TypeIn == 1 && ggc == 0
                [~, foundFrac1] = min(abs(pulseIny(VMinIndex{c}(s)- ...
                    round((500E-9)/dN):VMinIndex{c}(s)) - VMinFraction1));
                [~, foundFrac2] = min(abs(pulseIny(VMinIndex{c}(s)- ...
                    round((500E-9)/dN):VMinIndex{c}(s)) - VMinFraction2));
                
                foundFrac1 = foundFrac1 + VMinIndex{c}(s) + round(500E-9/dN);
                foundFrac2 = foundFrac2 + VMinIndex{c}(s) + round(500E-9/dN);
                
            else
                %}
            if scintPulse == 1
                'SCINT PULSE!!!!!!!'
                try
                    [~, foundFrac1] = min(abs(pulseIny(VMinIndex{c}(s)-round(100E-9/dN): ...
                        VMinIndex{c}(s)) - VMinFraction1));
                    [~, foundFrac2] = min(abs(pulseIny(VMinIndex{c}(s)-round(100E-9/dN): ...
                        VMinIndex{c}(s)) - VMinFraction2));
                    
                    foundFrac1 = foundFrac1 + VMinIndex{c}(s) - round(100E-9/dN);
                    foundFrac2 = foundFrac2 + VMinIndex{c}(s) - round(100E-9/dN);
                    
                catch
                    [~, foundFrac1] = min(abs(pulseIny(1: ...
                        VMinIndex{c}(s)) - VMinFraction1));
                    [~, foundFrac2] = min(abs(pulseIny(1: ...
                        VMinIndex{c}(s)) - VMinFraction2));
                    
                    %foundFrac1 = foundFrac1 + VMinIndex{c}(s) - round(10E-9/dN);
                    %foundFrac2 = foundFrac2 + VMinIndex{c}(s) - round(10E-9/dN);
                    
                end
                
            else
                
                earlyOne = smooth(pulseIny(VMinIndex{c}(s)-round(5E-6/dN): ...
                    VMinIndex{c}(s) + round(5E-6/dN)), 50);
                
                leap = diff(smooth(diff(earlyOne), 50));
                
                [~, risingEdge] = max(leap(200:end-200));
                [~, turningEdge] = min(leap(200:end-200));
                
                
                %plot(leap)
                
                
                %leap = diff(smooth(diff(smooth(pulseIny(VMinIndex{c}(s) - round(1500E-9/dN): ...
                %    VMinIndex{c}(s)), 100)), 100));
                
                %[~, risingEdge] = max(leap(200:end - 100));
                
                risingEdge = risingEdge + VMinIndex{c}(s) - round(5E-6/dN) + 200;
                turningEdge = turningEdge + VMinIndex{c}(s) - round(5E-6/dN) + 200;
                
                
                
                [~, midPulse] = min(abs(pulseIny(risingEdge: ...
                    turningEdge) - 0.20*VMin{c}(s)));
                
                midPulse = midPulse + risingEdge;
                
                
               
                %size(pulseIny)
                %size(leap)
                %plot(leap(200:end-100))
                
                %{
                [fracVal1, foundFrac1] = min(abs(pulseIny(VMinIndex{c}(s)-round(700E-9/dN): ...
                    VMinIndex{c}(s)) - VMinFraction1));
                [fracVal2, foundFrac2] = min(abs(pulseIny(VMinIndex{c}(s)-round(700E-9/dN): ...
                    VMinIndex{c}(s)) - VMinFraction2));
                
                foundFrac1 = foundFrac1 + VMinIndex{c}(s) - round(700E-9/dN);
                foundFrac2 = foundFrac2 + VMinIndex{c}(s) - round(700E-9/dN);
                %}
                
                [fracVal1, foundFrac1] = min(abs(pulseIny(risingEdge: ...
                    midPulse) - VMinFraction1));
                [fracVal2, foundFrac2] = min(abs(pulseIny(risingEdge: ...
                    midPulse) - VMinFraction2));
                
                foundFrac1 = foundFrac1 + risingEdge;
                foundFrac2 = foundFrac2 + risingEdge;
                
                %{
                 'ASDFIF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                risingEdge
                turningEdge
                midPulse
                VMinIndex{c}(s)
                foundFrac1
                foundFrac2
                pause(5)
                
                
                
                fracVal1
                VMinFraction1 + 0.25*VMinFraction1
                fracVal2
                fracVal1 
                        foundFrac2
                        foundFrac1
                %}
                if fracVal1 >= VMinFraction1 + 0.25*VMinFraction1 || fracVal2 < fracVal1 || ...
                        foundFrac2 < foundFrac1
                    'HERE'
                    [fracVal1, foundFrac1] = min(abs(pulseIny(VMinIndex{c}(s) - round(5000E-9/dN): ...
                        midPulse) - VMinFraction1));
                    [fracVal1, foundFrac2] = min(abs(pulseIny(VMinIndex{c}(s) - round(5000E-9/dN): ...
                        midPulse) - VMinFraction2));
                    
                    foundFrac1 = foundFrac1 + VMinIndex{c}(s) - round(5000E-9/dN);
                    foundFrac2 = foundFrac2 + VMinIndex{c}(s) - round(5000E-9/dN);
                    
                    %{
                    [~, foundFrac1] = min(abs(pulseIny(midPulse - round(100E-9/dN): ...
                        midPulse) - VMinFraction1));
                    [~, foundFrac2] = min(abs(pulseIny(midPulse - round(100E-9/dN): ...
                        midPulse) - VMinFraction2));
                    
                    foundFrac1 = foundFrac1 + VMinIndex{c}(s) - midPulse;
                    foundFrac2 = foundFrac2 + VMinIndex{c}(s) - midPulse;
                    %}
                        
                end
                
                
                if pulseInx(foundFrac2) - pulseInx(foundFrac1) > 500E-9 || ...
                        foundFrac2 < foundFrac1 || ...
                        fracVal1 >= VMinFraction1 + 0.25*VMinFraction1
                    
                    crossTime(s) = NaN;
                    riseTime = 1;
                    numPeaks = numPeaks - 1;
                    'Bad ELET'
                    continue;
                    
                end
                
            end
            
                %{
                if foundFrac1 == foundFrac2
                figure
                plot((pulseIny(VMinIndex{c}(s)-round(250E-9/dN):VMinIndex{c}(s))))
                pause(5)
                end
                  %}
                
                %{
                if c == 1
                    'foundFracs c = 1'
                    foundFrac1
                    foundFrac2
                    if foundFrac1 == foundFrac2
                        plot(pulseIny)
                        
                    end
                    
                end
                %}
                
                if findingNoise == true
                    foundNoise = false;
                    l = foundFrac2
                    
                    while foundNoise == false %|| foundFrac2 - l > round(100E-9/dN)
                        
                        if pulseIny(l-1) >= pulseIny(l)
                            foundFrac1 = l + 10;
                            VMinFraction1 = pulseIny(foundFrac1);
                            foundNoise = true;
                            
                        else
                            l = l - 1;
                            
                        end
                        
                    end
                    
                end
                
            %end
          
        catch
            crossTime(s) = NaN;
            riseTime = 1;
            numPeaks = numPeaks - 1;
            'ELET error'
            continue;
            
        end
        
    % ---------------------------------------------------------------------
    % Trapezoidal Shaping Timing
    % ---------------------------------------------------------------------
        
    elseif TypeIn == 2 && timingIn == 20 && T ~= 0
        %{
        stop1 = 0;
        stop2 = 0;
        
        while (stop1 == 0 || stop2 == 0) && VMinIndex{c}(s) > count2 && VMinIndex{c}(s) > count1 ...
                && count2 < 0.25*numel(pulseIny) && count1 < 0.25*numel(pulseIny) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% count2 was 3000
            % Going down peak edge by 'count' until less than VMinFraction.
            if stop1 == 0 && abs(pulseIny(VMinIndex{c}(s) - count1)) <= abs(VMinFraction1)
                foundFrac1 = (VMinIndex{c}(s) - count1);
                stop1 = 1;
            elseif stop1 == 0
                count1 = (count1 + 1);
            end
            
            % Going down peak edge by 'count' until less than VMinFraction.
            if  stop2 == 0 && abs(pulseIny(VMinIndex{c}(s) - count2)) <= abs(VMinFraction2)
                foundFrac2 = (VMinIndex{c}(s) - count2);
                stop2 = 1;
            elseif stop2 == 0
                count2 = (count2 + 1);
            end
        end
        
    % else if truth vector element is false (bad pulse).
    elseif T == 0
        noGood = true;
        riseTime = 1;
        crossTime(s) = 2;
        'bad truth 1'
        %}
        
        [~, foundFrac1] = min(abs(pulseIny(1:VMinIndex{c}(s)) - VMinFraction1));
        [~, foundFrac2] = min(abs(pulseIny(1:VMinIndex{c}(s)) - VMinFraction2));
        
    end
    
    if TypeIn == 3 && timingIn == 3 && T ~= 0
        try
            %{
            [output,~,maximum] = neuralNetworkFiltering(pulseIn,1000);
            
            [~,MF1Index] = min(abs(output-0.20*maximum));
            [~,MF2Index] = min(abs(output-0.75*maximum));
            
            x = linspace(MF2Index,MF1Index,MF1Index-MF2Index+1)';
            
            lastwarn('');
            tryFit = polyfit(x,output(MF2Index:MF1Index),1);
            [warnmsg,~] = lastwarn;
            %}
            vOutDiff1 = diff(pulseIny);
            vOutDiff1Smooth = smooth(vOutDiff1,251,'lowess');
            vOutDiff2 = diff(vOutDiff1Smooth);
            vOutDiff3Smooth = smooth(vOutDiff2,251,'lowess');
            %vOutDiff3 = diff(vOutDiff2Smooth);
            %vOutDiff3Smooth = smooth(vOutDiff3,251,'lowess');
            %figure
            %plot(vOutDiff3Smooth)
            
            trapezoidalMax = max(abs(vOutDiff3Smooth(100:numel(vOutDiff3Smooth)-1000)));
            
            halfTrapezoidalMax = trapezoidalMax/2;
            
            [~,index] = findpeaks(abs(vOutDiff3Smooth(1:numel(vOutDiff3Smooth)-1000)), ...
                'MinPeakHeight',halfTrapezoidalMax);
            
        catch
            warnmsg = 'bad';
            'warnmsg'
            
        end
    end
    
    if isempty(foundFrac) 
        foundFrac = 1;
        'foundFrac empty'
    end
    
    %{
    'TRUTHS'
    foundFrac1
    foundFrac2
    c
    (timingIn == 0 || timingIn == 2)
    foundFrac1 > round(0.5E-9/dN)
    foundFrac2 > round(0.5E-9/dN)
    count1 < 0.99*numel(pulseIny)
    count2 < 0.99*numel(pulseIny)
    foundFrac1 < 0.99*numel(pulseIny)
    foundFrac2 < 0.99*numel(pulseIny)
    foundFrac2 - foundFrac1 <= round(900E-9/dN)
    T ~= 0
    pulseIny(foundFrac2) < 0.99*VMin{c}(s)
    foundFrac1 ~= foundFrac2
    pulseIny(foundFrac2) - pulseIny(foundFrac1) >= 0.25*(VMinFraction2 - VMinFraction1) 
    %}
    
    % CFD pulseIn interpolation.
    %{
    'TRUTHS 2'
    c
    timingIn == 1 
    timingIn == 2 
    foundFrac > 1 
    count < 0.99*numel(pulseIny) 
    foundFrac < 0.99*numel(pulseIny) 
    T ~= 0 
    
    foundFrac
    %}
    if (TypeIn ~= 1 && s > 1) || (ggc == 1 && TypeIn == 1 && s > 1)
        'TypeIn & s 2'
        break;
        
    elseif (timingIn == 1 || timingIn == 2 || (timingIn == 4 && cfdANN == true)) && foundFrac > 1 ...
            && count < 0.99*numel(pulseIny) && foundFrac < 0.99*numel(pulseIny) && ...
            T ~= 0
        
        % Interpolating between values closest to the true fraction.
        
        % Creating vectors for interpolation between the
        % values closest to the true fraction.
        % -----------------------------------------------------------------------------------------------------
        'CFD foundFrac';
        time2 = pulseInx(foundFrac);
        
        if foundFrac > VMinIndex{c}(s)
            time1 = pulseInx(foundFrac+1);
            %voltage1 = pulseIny(foundFrac+1);
            interpChannel = linspace(foundFrac,foundFrac+1);
            interpTime = linspace(time2,time1);
        else
            time1 = pulseInx(foundFrac-1);
            %voltage1 = pulseIny(foundFrac-1);
            interpChannel = linspace(foundFrac-1,foundFrac);
            interpTime = linspace(time1,time2);
        end
        %voltage2 = pulseIny(foundFrac);
        
        % --------------------------------------------------------------------------------------------------
        interpolation = interp1(pulseIny,interpChannel,'linear');
        
        [~, fractionIndex] = min(abs(interpolation - VMinFraction));
        'CFD CrossTime';
        crossTime(s) = vpa(interpTime(fractionIndex));
        noGood = false;
        
    % --------------------------------------------------------------------------------------
    % Extrapolated leading edge timing discrimination for pulseIn.
    % --------------------------------------------------------------------------------------
        
    elseif (timingIn == 0 || timingIn == 2 || timingIn == 4) && foundFrac1 > round(0.5E-9/dN) ...
            && foundFrac2 > round(0.5E-9/dN) && ...
            count1 < 0.99*numel(pulseIny) && count2 < 0.99*numel(pulseIny) && ...
            foundFrac1 < 0.99*numel(pulseIny) && foundFrac2 < 0.99*numel(pulseIny) && ...
            pulseIny(foundFrac2) < 0.99*VMin{c}(s)
            %foundFrac2 - foundFrac1 <= round(1000E-9/dN) && T ~= 0 && ...
            % && foundFrac1 ~= foundFrac2 && ...
            %pulseIny(foundFrac2) - pulseIny(foundFrac1) >= 0 ...%0.1*(VMinFraction2 - VMinFraction1) 
            %&& (activeNeuron == 0 || exact(activeNeuron) == 1)
        
            'FOUNDFRAC1 AND 2';
            foundFrac1;
            foundFrac2;
            if foundFrac2 == foundFrac1
                foundFrac2 = foundFrac2 + 1;
                
            end
            
            %{
            if foundFrac2 - foundFrac1 > 900E-9/dN
                foundFrac2 = foundFrac1 + round(900E-9/dN);
                
            end
            %}
        % Determining the x (time) and y (voltage) values of the
        % discovered fractions.
        %fractionTime1 = pulseInx(foundFrac1);
        %fractionTime2 = pulseInx(foundFrac2);
        %fractionVolt1 = pulseIny(foundFrac1);
        %fractionVolt2 = pulseIny(foundFrac2);
        %{
        % Finding pulse edge to determine baseline.
        vOutDiff1 = diff(pulseIny);
        vOutDiff1Smooth = smooth(vOutDiff1,151,'moving');
        vOutDiff2 = diff(vOutDiff1Smooth);
        
        [~,timingIndex] = max(vOutDiff1);
        trapezoidalMax = max(abs(vOutDiff2(1:numel(vOutDiff2)-1)));
        halfTrapezoidalMax = trapezoidalMax/2;
        
        [~,index] = findpeaks(abs(vOutDiff2(10:numel(vOutDiff2)-10)), ...
            'MinPeakHeight',halfTrapezoidalMax);
        index = index + 10;
        %}
        try
            baseline = 0;%mean(pulseIny(round(13E-6/dN):round(14E-6/dN)));
        catch
            baseline = 0;
        end
        
        % Linear fitting to determine the x crossover value, which
        % is taken as the pick-off time for the gamma pulseIn.
        
        %{
        % -------------------------------------------------------------
        % Attempt at timing along fastest rising slope region of pulse.
        % -------------------------------------------------------------
        lastwarn('');
        pFit = polyfit(pulseInx(timingIndex-round(5E-9/dN): ...
            timingIndex+round(5E-9/dN)),pulseIny(timingIndex-round(5E-9/dN): ...
            timingIndex+round(5E-9/dN)),1);
        [warnMsg, ~] = lastwarn;
        %}
        
        %size(pulseInx(foundFrac1:foundFrac2))
        %size(pulseIny(foundFrac1:foundFrac2))
        try
            lastwarn('');
            pFit = polyfit(pulseInx(foundFrac1:foundFrac2),pulseIny(foundFrac1:foundFrac2),1);
            [warnMsg, ~] = lastwarn;
            
            if ~isempty(warnMsg)
                'WARNING'
                pFit = polyfit(pulseInx(foundFrac1:foundFrac2 + 1),pulseIny(foundFrac1:foundFrac2 + 1),1);
                
            end
            
            fitEval = polyval(pFit,pulseInx(foundFrac1:foundFrac2));
            baseline = fitEval(1);
            
            
            
            % Rise time filter based on linear fit.
            if TypeIn == 3
                
                riseTime = abs((VMinFraction2 - VMinFraction1)/pFit(1));
                
                goodRiseTime = true;
                %{
            if riseTime <= RTFUpIn && riseTime >= RTFLIn
                goodRiseTime = true;
                
            else
                goodRiseTime = false;
                'goodRiseTime false'
                
            end
                %}
            else
                goodRiseTime = true;
                
            end
        catch
            'fit error'
            
        end
        %valuation = polyval(pFit,pulseInx(foundFrac1:foundFrac2));
        
        % Ensuring no fit error when calculating crossTime.
        if goodRiseTime == true && ...
                (timingIn == 0 || timingIn == 2 || timingIn == 4)
            
            crossTime(s) = -(baseline + pFit(2))/pFit(1);
            noGood = false;
            
            %{
            if TypeIn == 1
            assignin('base','crossTime',crossTime(s));
            evalin('base','crossTimes = horzcat(crossTimes,crossTime);');
            evalin('base','flags = 1;');
            end
            %}
            % If truth vector (T) is not satisfied.
        else
            crossTime(s) = NaN;
            numPeaks = numPeaks - 1;
            'rise time 1'
            
        end
        
    elseif timingIn ~= 3
        noGood = true;
        crossTime(s) = NaN;
        riseTime = 1;
        'timingIn 4'
        c
        
    end
    
    if TypeIn == 3 && timingIn == 5 && T ~= 0 && 1 == 0% && isempty(output) == 0
        try
            if goodRiseTime == true% && isempty(warnmsg)
                %'dfsdsffds'
                %{
                numberOfNetworks = 1;
                networkOutput{1} = myNeuralNetworkFunction9(tryFit(1));
                %networkOutput{2} = myNeuralNetworkFunction4(output);
                %networkOutput{3} = myNeuralNetworkFunction8(output);
                %networkOutput{4} = myNeuralNetworkFunction8(output);
                %networkOutput{5} = myNeuralNetworkFunction9(output);
                
                networkOutputAverage = (networkOutput{1}(1))*maximum;% + ...
                    %networkOutput{2}(1) + networkOutput{3}(1))/3;% + ...
                    %networkOutput{3}(1) + networkOutput{4}(1) + ...
                    %networkOutput{5}(1))/numberOfNetworks;
                
                [~,crossing] = min(abs(pulseIn.y(100:VMinIndex{c}+1000) - ...
                    networkOutputAverage));
                %}
                %{
                time2 = pulseInx(index(1));
                if networkOutputAverage > crossing
                    time1 = pulseInx(crossing+1);
                    interpChannel = linspace(crossing,crossing+1);
                    interpTime = linspace(time2,time1);
                else
                    time1 = pulseInx(crossing-1);
                    interpChannel = linspace(crossing-1,crossing);
                    interpTime = linspace(time1,time2);
                end
                
                % --------------------------------------------------------------------------------------------------
                interpolation = interp1(pulseIny,interpChannel,'spline');
                
                [~, fractionIndex] = min(abs(interpolation - VMinFraction));
                %}
                crossTime(s) = pulseInx(index(1));%vpa(interpTime(fractionIndex));
                noGood = false;
                
            else
                crossTime(s) = NaN;
                numPeaks = numPeaks - 1;
                riseTime = 1;
                'rise time 2'
                
            end
            
        catch
            crossTime(s) = NaN;
            numPeaks = numPeaks - 1;
            riseTime = 1;
            'error'
            
        end
        
    end
    
    
    % Plot CFD and ELET markers on pulses with more than one peak
    doesIt = evalin('base','exist(''ax'',''var'');');
    doesIt2 = evalin('base','exist(''ax2'',''var'');');
    
    if doesIt == 1
        isIt = evalin('base','isvalid(ax);');
    else
        isIt = 0;
    end
    
    if doesIt2 == 1
        isIt2 = evalin('base','isvalid(ax2);');
    else
        isIt2 = 0;
    end
    
    if isIt == 0
        evalin('base','f = figure;');
        evalin('base','ax = subplot(1,1,1);');
        evalin('base','plotNext = false;');
        
    end
    
    if isIt2 == 0
        evalin('base','f2 = figure;');
        evalin('base','ax2 = subplot(1,1,1);');
        evalin('base','plotNext = false;');
        
    end
    
    if TypeIn == 1 && s == numPeaks && foundFrac ~= 0 && TypeIn == 1 && noGood == false
        assignin('base','pulseIny',pulseIny);
        assignin('base','foundFraction',foundFraction);
        %evalin('base','figure(ax)');
        evalin('base','plot(ax,pulseIny,''-p'',''MarkerIndices'',foundFraction,''MarkerFaceColor'',''red'',''MarkerSize'',15);');
        evalin('base','plotNext = true;');
        %{
        if numPeaks == 2
            pause(2)
        end
        %}
    end
    
    plotNext = evalin('base','plotNext;');
    size(pulseIny);
    size(pulseInx);
    
    if TypeIn == 3 && tempOn == false && timingIn == 0 && plotNext == true && foundFrac1 ~= 0 && foundFrac2 ~= 0 && TypeIn == 3 && noGood == false
        assignin('base','pulseIny',pulseIny);
        assignin('base','pulseInx',pulseInx);
        assignin('base','foundFrac1',foundFrac1);
        assignin('base','foundFrac2',foundFrac2);
        %evalin('base','figure(ax)');
        evalin('base','plot(ax2,pulseInx,pulseIny,''-p'',''MarkerIndices'',[foundFrac1 foundFrac2],''MarkerFaceColor'',''red'',''MarkerSize'',15);');
        evalin('base','plotNext = false;');
        
    end
    
    % ---------------------------------------------------------------------
    %{
    'FIRST'
    TypeIn == 3
    tempOn == true
    %}
    
    if TypeIn == 3 && tempOn == true %s == numPeaks && foundFrac ~= 0 && TypeIn == 1 && noGood == false
        assignin('base','pulseIny',pulseIny);
        assignin('base','foundFrac',foundFrac);
        %evalin('base','figure(ax)');
        evalin('base','plot(ax,pulseIny,''-p'',''MarkerIndices'', foundFrac,''MarkerFaceColor'',''red'',''MarkerSize'',15);');
        evalin('base','plotNext = true;');
        
        risingEdge
        VMinIndex{c}(s)
        
    end
    
    plotNext = evalin('base','plotNext;');
    %{
    'SECOND'
    TypeIn == 3 
    scintPulse == 1
    tempOn == false
    plotNext == true
    %}
    if TypeIn == 3 && tempOn == false && plotNext == true
        assignin('base','pulseIny',pulseIny);
        assignin('base','pulseInx',pulseInx);
        assignin('base','foundFrac',foundFrac);
        %evalin('base','figure(ax)');
        evalin('base','plot(ax2,pulseIny,''-p'',''MarkerIndices'', foundFrac,''MarkerFaceColor'',''red'',''MarkerSize'',15);');
        evalin('base','plotNext = false;');
        
    end
    %}
    %crossTime(s) = crossTime(s)/1E9;
end


% Breaking from the loop if there are more or fewer electrons detected
% than are set by the multi-stop conditions.
if TypeIn == 1 && ((multiStopIn == '1' && (numel(warnMsg) > 0) || ...
        (multiStopConditionsIn ~= 4 && multiStopConditionsIn ~= 5 ...
        && multiStopConditionsIn ~= numPeaks) ...
        || (multiStopConditionsIn == 5 && (numPeaks == 1 || numPeaks == 0)) ...
        || numPeaks > 5))
    
    noGood = true;
    riseTime = 1;
    numberOfPeaks = 0;
    
    %{
    if rejecting == true
    'NUMPEAKS'
        multiStopConditionsIn
        numPeaks
        pause(10)
    
    end
      %}  
    %{
        if integration < 0.005
            integration
            'integration 2'
            figure
            plot(pulseIny)
        end
    %}
    'numPeaks reject'
    
end

end