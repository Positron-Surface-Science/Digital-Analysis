function [amp,shapedPulse,noGood] = MWD(pulseIn,oldPulse,TFRFIn,TFTopIn)
format long
digits(15)

p = 1;
q = floor((1.25*p)/0.1);
Fs = (p/q)*pulseIn.desc.fs;
T = (1/Fs);
divid = round(288E-9/T);

L = round(TFRFIn/T)
G = round(TFTopIn/T)

amp = 0;

offset1 = -pulseIn.info.OFFSET;
offset = offset1 + 1100E-9;

is = 1;
plottingOn = 0;
plottingOn2 = 0;
plottingOn3 = 0;
plottingOn4 = 0;
amp = 0;
SLF = 0;
noGood = 0;
ANN = 0;

try
    
    measurePulse = smooth(pulseIn.y,101,'moving') - ...
        mean(pulseIn.y(round(offset/T-2E-6/T):round(offset/T-1E-6/T)));
    
    [m,ending] = max(measurePulse);
    [~,begin] = min(abs(measurePulse(round(-pulseIn.info.OFFSET/T-1E-6/T): ...
        round(-pulseIn.info.OFFSET/T))-0.1*m));
    
    begin = begin + round((-pulseIn.info.OFFSET/T-1E-6/T));
    
    riseTime = (pulseIn.x(round(ending*q/p)) - pulseIn.x(round(begin*q/p)));
    
    if riseTime > 1000E-9 || riseTime < 110E-9
        noGood = 1;
        amp = 0;
        shapedPulse = pulseIn.y;
        'rise time shaper'
        return;
        
    end

    if is == 1
        
        unNew = pulseIn.y;
        pulseDeconv.y = pulseIn.y;
        difference = pulseDeconv.y;
        shapedPulse = pulseDeconv.y;
        
        for n=2:numel(unNew)-1
            
            if n > round(offset/T)
                tau = 67545.59927*p/q;
                pulseDeconv.y(n) = pulseDeconv.y(n-1) + unNew(n) - unNew(n-1) + unNew(n-1)/tau;
                
            else
                pulseDeconv.y(n) = unNew(n);
            end
            
            
            if n > G + 10
                difference(n) = pulseDeconv.y(n) - pulseDeconv.y(n-G);
            else
                difference(n) = 0;
            end
            
            difference(n) = difference(n)*100;
            shapedPulse(n) = difference(n);
            
        end
        
        shapedPulse = shapedPulse(10:numel(shapedPulse)-10);
       
        indexBegin = round(offset/T) + round(L)%-round(100E-9/T);%round(15E-6/T)+round(G/2)% - round(1E-6/T);%index(2) + round(1000E-9/T);
        indexEnd = round(offset/T) + round(G-1200E-9/T)%round((G - round(1E-6/T))/divid)*divid;%round(15E-6/T)+round(G/2 + 10E-9/T)%index(3) - round(100E-9/T);
        
        if SLF == 1
            
            if RSS <= 0.0018 && p(1) <= 10E-6 && p(1) > -5.7E-5
                amp = values(maxIndex2)*100;
                shapedPulse = pulseDeconv.y;
                
            else
                amp = 0;
                shapedPulse = pulseIn.y;
                
            end
           
        end
        
        baseline = mean(pulseIn.y(round((offset1-999E-9)/T):round((offset1-500E-9)/T)));
        
        amp = mean(shapedPulse(indexBegin:indexEnd));
        
        if plottingOn == 1
           
            if amp >= 1.04736 && amp <= 1.08032
                'AMP'
                amp
                %figure
                %plot(pulseIn.y);
                %title('low side');
                assignin('base','p1',oldPulse(round(10E-6/T):round(35E-6/T)));
                evalin('base','p1s = horzcat(p1s,p1);');
                %evalin('base','baseline1 = horzcat(baseline1,baseline);');
                
            elseif amp >= 3.65295 && amp <= 3.67676
                'AMP'
                amp
                assignin('base','p2',oldPulse(round(10E-6/T):round(35E-6/T)));
                evalin('base','p2s =  horzcat(p2s,p2);');
                %evalin('base','baseline2 = horzcat(baseline1,baseline);');
                
            elseif amp >= 3.99902 && amp <= 4.03564
                'AMP'
                amp
                %assignin('base','pulse',pulseIn.y);
                %evalin('base','pulseHigh = pulseHigh + pulse;');
                %evalin('base','m=m+1');
                assignin('base','p3',oldPulse(round(10E-6/T):round(35E-6/T)));
                evalin('base','p3s =  horzcat(p3s,p3);');
                %evalin('base','baseline2 = horzcat(baseline1,baseline);');
                
            elseif amp > 4.70398 && amp <= 4.74609
                'AMP'
                amp
                assignin('base','p4',oldPulse(round(10E-6/T):round(35E-6/T)));
                evalin('base','p4s =  horzcat(p4s,p4);');
            
            elseif amp > 5.07935 && amp <= 5.10681
                'AMP'
                amp
                assignin('base','p5',oldPulse(round(10E-6/T):round(35E-6/T)));
                evalin('base','p5s =  horzcat(p5s,p5);');
                
            end
            
        end
        
        if plottingOn2 == 1
            whichBin = floor(amp/0.07);
            
            if whichBin > 0 && whichBin <= 100
                assignin('base','whichBin',whichBin);
                assignin('base','pulse',pulseIn.y(round(12E-6/T):round(24E-6/T)));
                %plot(shapedPulse(round(13E-6/T):round(22E-6/T)))
                evalin('base','geHistogram{whichBin} =  horzcat(geHistogram{whichBin},pulse);');
                
            end
            
        end
        
        if plottingOn3 == 1
            
            bin3Size = (7/4096);
            binSize = 0.0001;
            
            whichBinGe = floor(amp/bin3Size);
            whichBinBL = floor((baseline+0.004)/binSize);
            assignin('base','whichBinGe',whichBinGe);
            assignin('base','whichBinBL',whichBinBL);
            
            if whichBinGe > 0 && whichBinBL > 0 && whichBinBL <= 200 ...
                    && whichBinGe <= 4096
                
                evalin('base','geVsBL(whichBinGe,whichBinBL) = geVsBL(whichBinGe,whichBinBL) + 1;');
                
            end
            
        end
        
        if plottingOn4 == 1
            
            plot(pulseIn.y(round(10E-6/T):round(35E-6/T)))
            assignin('base','trace',pulseIn.y(round(10E-6/T):round(35E-6/T)));
            evalin('base','traces =  horzcat(traces,trace);');
            
        end
        
    elseif is == 2
        for n=2:numel(pulseIn.y)-1
            pulseDeconv.y(n) = pulseDeconv.y(n-1) + pulseIn.y(n) - pulseIn.y(n-1) + pulseIn.y(n-1)/tau;
            
            if n > 1
                difference(n) = pulseDeconv.y(n) - pulseDeconv.y(n-1);
            else
                difference(n) = 0;
            end
            
        end
        
        shapedPulse = smooth(difference(100:numel(difference)-500),G,'moving');
        shapedPulse = shapedPulse(100:numel(shapedPulse)-500)*1E6;
        
        vOutDiff1 = diff(shapedPulse);
        vOutDiff1Smooth = smooth(vOutDiff1,501,'lowess');
        vOutDiff2 = diff(vOutDiff1Smooth);
        
        trapezoidalMax = max(abs(vOutDiff2(1000:numel(vOutDiff2)-1000)));
        
        halfTrapezoidalMax = trapezoidalMax/2;
        %plot(vOutDiff2)
        [~,index] = findpeaks(abs(vOutDiff2(1000:numel(vOutDiff2)-1000)), ...
            'MinPeakHeight',halfTrapezoidalMax);
        
        index = index + 1000;
        
        indexBegin = index(2);
        indexEnd = index(3);
        
        amp = mean(shapedPulse(indexBegin:indexEnd)) - mean(shapedPulse(1000:indexBegin-100));
        
    elseif is == 3
        
         for n=2:numel(pulseIn.y)-1
            
            if n > G
                difference(n) = pulseIn.y(n) - pulseIn.y(n-G) + (1/tau)*sum(pulseIn.y(n-G:n-1));
            else
                difference(n) = 0;
            end
            
            if n > L
                trap(n) = (1/L)*sum(difference(n-L:n-1));
            end
            
         end
        
        vOutDiff1 = diff(pulseIn.y);
        vOutDiff1Smooth = smooth(vOutDiff1,101,'lowess');
        vOutDiff2 = diff(vOutDiff1Smooth);
        
        trapezoidalMax = max(abs(vOutDiff2(100:numel(vOutDiff2)-100)));
        
        halfTrapezoidalMax = trapezoidalMax/2;
        %plot(vOutDiff2)
        [~,index] = findpeaks(abs(vOutDiff2(100:numel(vOutDiff2)-100)), ...
            'MinPeakHeight',halfTrapezoidalMax);
        
        index = index + 100;
        
        index = index - round(50E-9/pulseIn.desc.Ts);
        plot(pulseIn.y(12000:index(1)+1000),'-p','MarkerIndices',index(1)-12000, ...
            'MarkerFaceColor','red','MarkerSize',15);
        shapedPulse = trap(100:numel(difference)-100)*1E3;
        amp = max(shapedPulse) - mean(pulseIn.y(100:10000));%sum(shapedPulse(index(1):index(2)));%+round((100E-9/pulseIn.desc.Ts))));
        
    elseif is == 4
        undo = log(pulseIn.y);
        maximum3 = max(undo);
        assignin('base','maximum3',maximum3);
        assignin('base','undo',undo);
        
        unNew = evalin('base','undo(20000:numel(undo),1) - (maximum3/maximum2)*unnamed(3001:numel(unnamed),1);');
        unNew = vertcat(pulseIn.y(1:19999),real(exp(unNew)));
        pulseDeconv = unNew;
        
        for n=2:numel(unNew)-1
            tau = 1/(3.44699E-5);%;/(T/(10.^(-9)));
            pulseDeconv(n) = pulseDeconv(n-1) + unNew(n) - unNew(n-1) + unNew(n-1)/tau;
            
        end
        
        shapedPulse = pulseDeconv*100;
        %plot(shapedPulse)
        %mean(pulseIn.y(11000:14000))*100
        %mean(shapedPulse(10000:30000))
        if min(shapedPulse(1:numel(shapedPulse)-10)) > 0
            amp = mean(shapedPulse()) - mean(shapedPulse(11000:14000))*100;
        else
            amp = 0;
        end
        
    elseif is == 5
        assignin('base','pulse',pulseIn.y);
        evalin('base','pulses = pulses + pulse;');
        amp = 0;
        shapedPulse = [];
        
    end
    
catch
    shapedPulse = pulseIn.y;
    noGood = 1;
    'ERROR'
    
end

end
%}