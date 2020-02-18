function [crossTime,noGood,numberOfPeaks,riseTime] = paes(timingIn,TypeIn,c,...
    multiStopIn,multiStopConditionsIn,numPeaks,pulseIn,RTFLIn,RTFUpIn,ELETMatrix,...
    DiscriminationVector,neuron,activeNeuron,ANN,ggc)

% --------------------------------------------------------------------------
% 'PAES' function which calculates the crossing "pick-off" time of detector
% pulses.
% --------------------------------------------------------------------------

pulseIny = pulseIn.y;
pulseInx = pulseIn.x;
dN = pulseIn.desc.Ts; % Time per channel
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
a = fieldnames(DiscriminationVector);
exact = getfield(DiscriminationVector,a{1});
% ----------------------------------
exact(exact == 0) = 1;
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
        'MinPeakDistance',round((5*10.^(-9))/dN), ...
        'MinPeakHeight',1*10.^(-3),'MinPeakProminence',0.8*10.^(-3), ...
        'MinPeakWidth',round((5*10.^(-9))/dN),'WidthReference', ...
        'halfheight','MaxPeakWidth',round((40*10.^(-9)/dN)));
    
    [warnMsg, ~] = lastwarn;
    VMinIndex{c} = VMinIndex{c} + 50;
    numPeaks = numel(VMinIndex{c});
    numberOfPeaks = numPeaks;
    crossTime = NaN(numPeaks,1);
    'numPeaks'
    numPeaks
    
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
        multiStopConditionsIn = 3;
        
    elseif multiStopIn == '1' && multiStopConditionsIn == 7 && numPeaks == 3
        multiStopConditionsIn = 3;
        
    elseif multiStopIn == '1' && multiStopConditionsIn == 6
        multiStopConditionsIn = 2
        
    end
   
    if multiStopIn == '0'
        numPeaks = 1;
        numberOfPeaks = 1;
        
    end
    
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
    
end

noGood = false;
VMinFraction = zeros(1,numPeaks);
VMinFraction1 = zeros(1,numPeaks);
VMinFraction2 = zeros(1,numPeaks);

numberOfPeaks = numPeaks;

for s=1:numberOfPeaks
    
    if isempty(activeNeuron) ||  ...
            isempty(VMinIndex{c}) || ...
            VMinIndex{c}(1) - 100 <= 0
        crossTime(s) = NaN;
        numPeaks = 0;
        'neuron error 1'
        break;
    
    end
        
    if TypeIn == 3 && ANN == 1 && activeNeuron == 0
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
    
    % Integrating peak and region before peak to reduce periodic noise
    % pulses.
    if (TypeIn == 1 && ggc == 0) && VMinIndex{c}(s)-round((100*10.^(-9)/dN)) > 0
        integration = trapz(pulseIny(VMinIndex{c}(s)-round((100*10.^(-9)/dN)) ...
            :VMinIndex{c}(s)));
    elseif TypeIn == 1
        integration = 0;
        %'TypeIn 1'
    end
    
    if (TypeIn == 1 && ggc == 0) && integration < 0.005
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
        VMinFraction = (0.42*VMin{c}(s));
        VMinFraction1 = (0.07*VMin{c}(s));
        VMinFraction2 = (0.09*VMin{c}(s));
        
        if ggc == 1 && TypeIn == 1
            VMinFraction1 = (0.10*VMin{c}(s));
            VMinFraction2 = (0.30*VMin{c}(s));
            
        elseif TypeIn == 3 && ANN == 1
            VMinFraction1 = (ELETMatrix(activeNeuron,1)*0.01*VMin{c}(s));
            VMinFraction2 = (ELETMatrix(activeNeuron,2)*0.01*VMin{c}(s));
            VMinFraction = (((ELETMatrix(activeNeuron,1)+ELETMatrix(activeNeuron,2)))/2)*0.01*VMin{c}(s);
            'ELET MATRIX 1'
            ELETMatrix(activeNeuron,1);
            ELETMatrix(activeNeuron,2);
            %VMinFraction1 = (0.07*VMin{c}(s));
            %VMinFraction2 = (0.09*VMin{c}(s));
            
        elseif ggc == 0 && TypeIn == 1
            VMinFraction1 = (0.1*VMin{c}(s));
            VMinFraction2 = (0.3*VMin{c}(s));
            
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
    
    elseif  (timingIn == 1 || timingIn == 2) && T ~= 0 && VMinIndex{c}(s) > round(700E-9/dN)
        
        
        if TypeIn == 1 && ggc == 0
            [~, foundFrac] = min(abs(pulseIny(VMinIndex{c}(s)- ...
                round(50*10.^(-9)/dN):VMinIndex{c}(s)) - VMinFraction));
            
            foundFrac = foundFrac + (VMinIndex{c}(s) - round(50*10.^(-9)/dN));
            
            
        else
            [~, foundFrac] = min(abs(pulseIny(VMinIndex{c}(s)-round(700E-9/dN): ...
                    VMinIndex{c}(s)) - VMinFraction));
            foundFrac = foundFrac + (VMinIndex{c}(s) - round(700E-9/dN));
            
        end
        
    % ---------------------------------
    % ELET: finding channels corresponding to the two fractions.
    % ---------------------------------
    
    elseif (timingIn == 0) && T ~= 0 && ...
            VMinIndex{c}(s) > round(900E-9/dN)
            
        
        try
            
                [~, foundFrac1] = min(abs(pulseIny(VMinIndex{c}(s)-round(900E-9/dN): ...
                    VMinIndex{c}(s)) - VMinFraction1));
                [~, foundFrac2] = min(abs(pulseIny(VMinIndex{c}(s)-round(900E-9/dN): ...
                    VMinIndex{c}(s)) - VMinFraction2));
                
                foundFrac1 = foundFrac1 + VMinIndex{c}(s) - round(900E-9/dN);
                foundFrac2 = foundFrac2 + VMinIndex{c}(s) - round(900E-9/dN);
                
                if findingNoise == true
                    foundNoise = false;
                    l = foundFrac2
                    
                    while foundNoise == false
                        
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
            vOutDiff1 = diff(pulseIny);
            vOutDiff1Smooth = smooth(vOutDiff1,251,'lowess');
            vOutDiff2 = diff(vOutDiff1Smooth);
            vOutDiff3Smooth = smooth(vOutDiff2,251,'lowess');
            
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
    
    % CFD pulseIn interpolation.
    
    if (TypeIn ~= 1 && s > 1) || (ggc == 1 && TypeIn == 1 && s > 1)
        'TypeIn & s 2'
        break;
        
    elseif (timingIn == 1 || timingIn == 2) && foundFrac > round(5*10.^(-9)/dN) ...
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
        
        % --------------------------------------------------------------------------------------------------
        interpolation = interp1(pulseIny,interpChannel,'linear');
        
        [~, fractionIndex] = min(abs(interpolation - VMinFraction));
        'CFD CrossTime';
        crossTime(s) = vpa(interpTime(fractionIndex));
        noGood = false;
        
    % --------------------------------------------------------------------------------------
    % Extrapolated leading edge timing discrimination for pulseIn.
    % --------------------------------------------------------------------------------------
        
    elseif (timingIn == 0 || timingIn == 2) && foundFrac1 > round(5E-9/dN) ...
            && foundFrac2 > round(5E-9/dN) && ...
            count1 < 0.99*numel(pulseIny) && count2 < 0.99*numel(pulseIny) && ...
            foundFrac1 < 0.99*numel(pulseIny) && foundFrac2 < 0.99*numel(pulseIny) && ...
            pulseIny(foundFrac2) < 0.99*VMin{c}(s)
        
            if foundFrac2 == foundFrac1
                foundFrac2 = foundFrac2 + 1;
                
            end
            
        try
            baseline = 0;
        catch
            baseline = 0;
        end
        
        % Linear fitting to determine the x crossover value, which
        % is taken as the pick-off time for the gamma pulseIn.
        
        % -------------------------------------------------------------
        % Attempt at timing along fastest rising slope region of pulse.
        % -------------------------------------------------------------
        
        try
        lastwarn('');
        pFit = polyfit(pulseInx(foundFrac1:foundFrac2),pulseIny(foundFrac1:foundFrac2),1);
        [warnMsg, ~] = lastwarn;
        
        
            fitEval = polyval(pFit,pulseInx(foundFrac1:foundFrac2));
            baseline = fitEval(1);
            
        
        
        % Rise time filter based on linear fit.
        if TypeIn == 3
            
            riseTime = abs((VMinFraction2 - VMinFraction1)/pFit(1));
            
            if riseTime < RTFUpIn && riseTime > RTFLIn
                goodRiseTime = true;
                
            else
                goodRiseTime = false;
                'goodRiseTime false'
                
            end
            
        else
            goodRiseTime = true;
            
        end
        catch
            'fit error'
            
        end
        
        % Ensuring no fit error when calculating crossTime.
        if goodRiseTime == true && isempty(warnMsg) && ...
                (timingIn == 0 || timingIn == 2)
            
            crossTime(s) = -(baseline + pFit(2))/pFit(1);
            noGood = false;
            
        else
            crossTime(s) = NaN;
            numPeaks = numPeaks - 1;
            'rise time 1'
            
        end
        
    elseif timingIn ~= 4
        noGood = true;
        crossTime(s) = NaN;
        riseTime = 1;
        'timingIn 4'
        c
        
    end
    
    if TypeIn == 3 && timingIn == 4 && T ~= 0 && 1 == 0
        try
            if goodRiseTime == true
               
                crossTime(s) = pulseInx(index(1));
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
    
    'numPeaks reject'
    
end

end