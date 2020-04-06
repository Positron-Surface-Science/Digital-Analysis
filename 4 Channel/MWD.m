function [amp,shapedPulse,noGood] = MWD(pulseIn,oldPulse,TFRFIn,TFTopIn)
format long
digits(15)

p = 1;%round(q/(L*pulseIn.desc.fs))
q = floor((1.25*p)/0.1);
Fs = (p/q)*pulseIn.desc.fs;
T = (1/Fs);%q/p)*pulseIn.desc.Ts
TFRFIn
divid = round(288E-9/T);%243E-9;
%gRound = round(TFTopIn/T);

L = round(TFRFIn/T) %3*10^(-6)
%G = round(gRound/divid)*divid
G = round(TFTopIn/T)

%{
best = evalin('base','best;');
activeNeuron = evalin('base','activeNeuron');
LT = L
GT = G
L = round((best(1,activeNeuron)*1E-6)/T)
G = round((best(2,activeNeuron)*1E-6)/T)
%}

if L == 0 && G == 0
    L = LT;
    G = GT;
    
end

amp = 0;

%{
for n=1:numel(pulseIn.y)-1
    assignin('base','pulseIn',pulseIn.y);
    assignin('base','n',n);
    ap = ap + evalin('base','h(n)*pulseIn(n);');
    
end

assignin('base','ap',ap);
pulseIn.y = real(evalin('base','ap.*s;'));
%}

offset1 = -pulseIn.info.OFFSET;
offset = offset1 + 1100E-9;

%plot(pulseIn.y)
%pulseIn.y = resample(pulseIn.y,p,q);

%tau = (58000/(T/(10.^(-9))));
%pulseDeconv.y = pulseIn.y;

%difference = pulseDeconv.y;
%shapedPulse = pulseDeconv.y;

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
    
%{
    if ANN == 1
        ap = 0;
        newvIn = pulseIn.y(round(12E-6/T):round(24E-6/T));
        for n=1:numel(newvIn)-1
            assignin('base','pulse',newvIn);
            assignin('base','n',n);
            ap = ap + evalin('base','h(n)*pulse(n);');
            
        end
        
        assignin('base','ap',ap);
        newestvIn = real(evalin('base','ap.*s;'));
        
        pulseIn.y = [zeros(1,round(12E-6/T)) newestvIn];
        
    end
    %}
    %[~,indexMax] = max(pulseIn.y);
    
    %pulseIn.y = pulseIn.y - mean(pulseIn.y(round(13E-6/T):round(14.5E-6/T)));
    
    if is == 1
        % -----------------------------------------------------------------
        % DNL Correction
        % -----------------------------------------------------------------
        %{
        pulse = pulseIn.y(round(11E-6/T):round(27E-6/T))';
        pulse2 = pulse;
        
        r = evalin('base','r');
        
        undone = log(pulse);
        r = r*max(real(undone));
        
        undone(550:numel(undone)) = undone(550:numel(undone)) - r';
        
        pulseIn.y = vertcat(pulseIn.y(1:round(11E-6/T)-1),real(exp(undone))',pulseIn.y(round(27E-6/T)+1:numel(pulseIn.y)));
        %}
        %{
        undone = log(pulseIn.y);
        %[~,maxIndex] = max(pulseIn.y);
        maxIndex = round(offset/T);%maxIndex + round(0.25E-6/T);
        endLinear = round(offset/T + 0.5E-6/T);
        %maxIndex = round(15.5E-6/T);%maxIndex + round(1E-6/T);
        %maximum3 = max(undo);
        %assignin('base','maximum3',maximum3);
        %assignin('base','undo',undo);
        %maxIndex+round(20E-6/T));%
        undo2 = undone(maxIndex:maxIndex+endLinear-1);
        %undo3 = undone(maxIndex+endLinear+1:numel(undone)-1);
        x = linspace(maxIndex-900E-9/T,maxIndex-900E-9/T+endLinear,endLinear)';
        %x7 = linspace(maxIndex+endLinear+1,numel(undone),numel(undone)-(maxIndex+endLinear+1))';
        %assignin('base','x',x)
        %undone = undo2;
        %unNew = evalin('base','undo(maxIndex:numel(undo),1) - (maximum3/maximum2)*unnamed(8001:numel(unnamed),1);');
        
        % -----------------------------------------------------------------
        % Fitting to discover tau for each pulse.
        % Linear fit produces better resolution than polynomial fit.
        % -----------------------------------------------------------------
        
        %undo = smooth(undo,501,'lowess');
        p = polyfit(x,real(undo2),3);
        %p7 = polyfit(x7,real(undo3),5);
        pV = polyval(p,x);
        
        amp = exp(pV(1))*100
        shapedPulse = undo2;
        noGood = 0;
        return;
        
        %plot(x,undo(25000:numel(undo)),x,unNew)
        %figure
        %plot(x,undo,x,polyval(p,x))
        %unNew = pulseIn.y;%vertcat(pulseIn.y(1:24999),real(exp(unNew)));%padarray(real(exp(unNew)),[5000 1],0,'pre');
        %pulseDeconv.y = unNew(1:numel(unNew));
        %}
        %pulseIn.y = pulseIn.y - mean(pulseIn.y(round(13.5E-6/T):round(14E-6/T)));
        %oldPulse = pulseIn.y;
        %baseline = mean(pulseIn.y(round((offset-999E-9)/T):round((offset-500E-9)/T)));
        %pulseIn.y = pulseIn.y - baseline;
        %assignin('base','pulseIn',pulseIn.y);
        unNew = pulseIn.y;
        pulseDeconv.y = pulseIn.y;
        difference = pulseDeconv.y;
        shapedPulse = pulseDeconv.y;
        
        for n=2:numel(unNew)-1%-1/(4*p(1)*n^3 + 3*p(2)*n^2 + 2*p(3)*n + p(4));
            %-1/(3*p(1)*n^2 + 2*p(2)*n + p(3));
            %{
            if n < 7235 && n >= 7235
                tau = -1/(9*7.66804E-27*n^8 - 3.43752E-22*n^7 + 7*5.80467E-18*n^6 - 6*3.37467E-14*n^5 - 5*2.72652E-10*n^4 + 4*6.19746E-6*n^3 - 3*0.04847*n^2 + 2*202.60999*n - 452187.83174);
            elseif n >= 7235
            %}
                %tau = -1/(-5*2.86981E-22*n^4 + 4*2.19799E-17*n^3 - 3*6.65158E-13*n^2 + 2*9.90794E-9*n - 1.07269E-4);
            %else
            
                %tau = (54000/T)/1E9;
            %end
            
            %tau = -1/p(1);
            
            if n > round(offset/T)%&& n <= maxIndex + endLinear%-1/p(1); %65145.59927*p/q;%
                tau = 67545.59927*p/q;%-1/p(1);%-1/(5*p(1)*n^4 + 4*p(2)*n^3 + 3*p(3)*n^2 + 2*p(4)*n + p(5));%-1/(4*p(1)*n^3 + 3*p(2)*n^2 + 2*p(3)*n + p(4));%-1/(3*p(1)*n^2 + 2*p(2)*n + p(3));%-1/(7*p(1)*n^6 + 6*p(2)*n^5 + 5*p(3)*n^4 + 4*p(4)*n^3 + 3*p(5)*n^2 + 2*p(6)*n + p(7));
                pulseDeconv.y(n) = pulseDeconv.y(n-1) + unNew(n) - unNew(n-1) + unNew(n-1)/tau;
                
            %elseif n > maxIndex + endLinear
            %    tau = -1/(5*p7(1)*n^4 + 4*p7(2)*n^3 + 3*p7(3)*n^2 + 2*p7(4)*n + p7(5));
            %    pulseDeconv.y(n) = pulseDeconv.y(n-1) + unNew(n) - unNew(n-1) + unNew(n-1)/tau;
            
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
        
        %plot(pulseDeconv.y)
        
        %{
        % ------------------------------------------
        undo = log(difference);
        [~,maxIndex] = max(difference);
        maxIndex = round((offset)/T);
        %maxIndex = maxIndex + round(0.25E-6/T);
        undo2 = undo(maxIndex:maxIndex+G-round(1E-6/T)-1);%numel(undo)-1);
        x = linspace(maxIndex,maxIndex+G-round(1E-6/T),G-round(1E-6/T))';
        undo = undo2;
        size(x)
        size(undo2)
        p = polyfit(x,real(undo),1);
        f = polyval(p,x);
        unNew2 = difference;
        pulseDeconv.y = unNew2(1:numel(unNew));
        for n=2:difference-1
            if n >= maxIndex
                tau = -1/p(1); %-1/(5*p(1)*n^4 + 4*p(2)*n^3 + 3*p(3)*n^2 + 2*p(4)*n + p(5));
                pulseDeconv.y(n) = pulseDeconv.y(n-1) + unNew2(n) - unNew2(n-1) + unNew2(n-1)/tau;
            else
                pulseDeconv.y(n) = unNew2(n);
            end
        end
        maxIndex
        round(offset/T)
        plot([f real(undo)])
        plot(pulseDeconv.y)
        % -------------------------------------------
        %}
        % -----------------------------------------------------------------
        % END
        % -----------------------------------------------------------------
            %{
                tau = -1/(-7.97737E-5 + 2*2.65368E-9*n - 3*7.63298E-14*n^2 + 4*1.07774E-18*n^3 - 5*5.96708E-24*n^4);
            %}
        %plot(shapedPulse)
        shapedPulse = shapedPulse(10:numel(shapedPulse)-10);
        %shapedPulse2 = shapedPulse;
        
        %shapedPulse = smooth(shapedPulse,L,'moving');
        %{
        vOutDiff1 = diff(shapedPulse);
        vOutDiff1Smooth = smooth(vOutDiff1,251,'lowess');
        vOutDiff2 = diff(vOutDiff1Smooth);
        
        trapezoidalMax = max(abs(vOutDiff2(100:numel(vOutDiff2)-1000)));
        
        halfTrapezoidalMax = trapezoidalMax/2;
        
        [~,index] = findpeaks(abs(vOutDiff2(1:numel(vOutDiff2)-1000)), ...
            'MinPeakHeight',halfTrapezoidalMax,'MinPeakDistance',0.005*(G));
        %}
        %plot(vOutDiff2)
        
        %plot(shapedPulse)
        %[~,in] = max(shapedPulse);
        
        indexBegin = round(offset/T) + round(L)%-round(100E-9/T);%round(15E-6/T)+round(G/2)% - round(1E-6/T);%index(2) + round(1000E-9/T);
        indexEnd = round(offset/T) + round(G-1200E-9/T)%round((G - round(1E-6/T))/divid)*divid;%round(15E-6/T)+round(G/2 + 10E-9/T)%index(3) - round(100E-9/T);
        
        %iB = round(offset/T) + round(LT);
        %iE = round(offset/T) + round(GT-1200E-9/T);
        
        %BLIndexStart = index(1) - round(500E-9/T);
        %BLIndexEnd = index(4) + round(500E-9/T);
        
        % -----------------------------------------------------------------
        % Testing baseline fitting
        % -----------------------------------------------------------------
        %BLBegin = mean(pulseDeconv.y(BLIndexStart-round(1000E-9/T):BLIndexStart));
        %BLEnd = mean(pulseDeconv.y(BLIndexStart-round(100E-9/T):BLIndexStart));
        
        %plot(x3,values2)
        % Straight line fit to deconvoluted pulse.
        %[~,maxIndex2] = max(pulseIn.y);
        %{
        maxIndex3 = round(offset/T);
        x1 = linspace(1,maxIndex3+G-round(1E-6/T),maxIndex3+G-round(1E-6/T))';
        %maxIndex2;% + round(100E-9/T);
        pulseFit1 = shapedPulse;
        %size(x1(maxIndex3:numel(x1)));
        %size(pulseFit1(maxIndex3:maxIndex3+G-round(1E-6/T)));
        
        [p2,S,mu] = polyfit(x1(maxIndex3:numel(x1)),pulseFit1(maxIndex3:maxIndex3+G-round(1E-6/T)),1);
        [values,delta] = polyval(p2,x1,S,mu);
        RSS = sum(delta.^2);
        %}
        if SLF == 1
            
            % -----------------------------------------------------------------
            % Testing straight line fit output
            % -----------------------------------------------------------------
            %[p3,S2,mu2] = polyfit(x2(1:BLIndexStart),pulseFit(1:BLIndexStart),1);
            %[values2,delta2] = polyval(p3,x2,S2,mu2);
            %plot([values values2 pulseFit])
            %assignin('base','values',values);
            %assignin('base','values2',values2);
            
            %values = values;% - values2;
            
            p(1);
            RSS;
            if RSS <= 0.0018 && p(1) <= 10E-6 && p(1) > -5.7E-5
                amp = values(maxIndex2)*100;
                shapedPulse = pulseDeconv.y;
                
            else
                amp = 0;
                shapedPulse = pulseIn.y;
                
            end
           
        end
        %{
            x2 = linspace(1,numel(shapedPulse)-10,numel(shapedPulse)-10)';
            pulseFit = shapedPulse(1:numel(shapedPulse)-10);
            [p3,S2,mu2] = polyfit(x2(G:round((offset-1E-6)/T)),pulseFit(G:round((offset-1E-6)/T)),1);
            [values2,delta2] = polyval(p3,x2,S2,mu2);
            
            RSS2 = sum(delta2.^2);
            %shapedPulse = shapedPulse(1:numel(shapedPulse)-10) - values2;
            %plot([values2 pulseFit])
        %}
        %'index'
        %numel(index)
        %baselineStartIndex = G;
        %((G-L))
        %plot(pulseDeconv.y);
        %amp = mean(shapedPulse(indexEnd:indexEnd+25)) - ...
        %mean(shapedPulse(baselineStartIndex:index(1)-10));
        %plot(shapedPulse);
        %abs(difference(indexEnd) - difference(indexBegin))/(indexEnd - indexBegin)
        %plot(shapedPulse(indexBegin:indexEnd));
        %figure
        %plot(vOutDiff2,'-p','MarkerIndices',[index(1) index(2) index(3) index(4)])
        
        %shapedPulse = smooth(shapedPulse,L,'moving');
        
        %abs(difference(indexEnd) - difference(indexBegin))/(indexEnd - indexBegin)
        %{
        RSS
        RSS2
        p2(1)
        p3(1)
        %}
        %shapedPulse = smooth(shapedPulse,G-L,'moving');
        
        %for l=1:5
        %    shapedPulse = smooth(shapedPulse,G-L,'moving');
        
        %end
        %{
        a = evalin('base','a');
        baseline = mean(mean(shapedPulse(round(12E-6/T):round(14E-6/T))) + mean(a));
        
        for t=2:5
            a(t-1) = a(t);
        end
        
        a(5) = baseline;
        assignin('base','a',a);
        %}
        
        baseline = mean(pulseIn.y(round((offset1-999E-9)/T):round((offset1-500E-9)/T)));
        %{
        x = baseline;
        y0 = 0.00661;
        A1 = -1.48569E-17;
        t1 = -0.00221;
        A2 = -8.67607E-4;
        t2 = -0.04968;
        A3 = -9.1457E-4;
        t3 = -0.05074;
        
        y = A1*exp(-x/t1) + A2*exp(-x/t2) + A3*exp(-x/t3) + y0 - 0.001;
         %}
        %indexBegin
        
        amp = mean(shapedPulse(indexBegin:indexEnd))
        %{
        amp2 = mean(shapedPulse(iB:iE))
        
        if amp2 >= 5.75 && amp >= 5.75 && amp2 ~= amp && best(2,activeNeuron) ~= 0
            amp = amp;% - (best(3,activeNeuron) - 6.58);
            
        elseif amp2 <= 5.75
            amp = 0;%amp - (6.6 - best(3,activeNeuron));%amp2;
            
        else
            amp = 0;
            
        end
        %}
        %{
        try
            if amp >= 6.2 && amp <= 7.3
            
                assignin('base','amp',amp);
                evalin('base','index = find(allNeurons.mwdNeurons{activeNeuron}(:,numberRuns) == 0);');
                evalin('base','allNeurons.mwdNeurons{activeNeuron}(index(1),numberRuns) = amp;');
                
            end
            
        catch
            'vector full'
            
        end
        %}
        %amp = amp - 0.5*amp*baseline^2;
        %assignin('base','pulseInx',pulseIn.x);
        %assignin('base','pulseIn',pulseIn.y);
        %assignin('base','pulseDeconv',pulseDeconv.y);
        %assignin('base','shapedPulse',shapedPulse);
        %amp = amp + (baseline*100)^2
        
        %beforeBaseline = mean(shapedPulse(round((offset-2.5E-6)/T):round((offset-1E-6)/T)));
        %afterBaseline = mean(shapedPulse(round((offset+1E-6)/T)+G:round((offset+2E-6)/T)+G));
        %{
        if baseline >= -0.00161 && baseline <= -0.00159
            amp = mean(shapedPulse(indexBegin:indexEnd));% - beforeBaseline);
            
        else
            amp = 0;
            
        end
        %}
        %plot([values(1:numel(pulseFit1)),pulseFit1])
        %plot(pulseIn.y)
        %numel(index)
        %abs(difference(indexEnd) - difference(indexBegin))/(indexEnd - indexBegin) < 2.5E-4 && ...
        %smoothed = smooth(diff(shapedPulse),101,'moving');
        %timeAboveDiff = sum(abs(pulseIn.x(smoothed(:) >= 3E-3)))
        %plot(shapedPulse(round((offset+1E-6)/T)+G:round((offset+2E-6)/T)+G))
        %plot(shapedPulse(indexBegin:indexEnd))
        %timeAboveDiff = 0.1;
        %if timeAboveDiff < 1
            %{
                RSS2 <= 7 && RSS <= 4 && RSS > 1E-5 && ... 
                abs(p2(1)) <= 7E-3 && abs(p3(1)) <= 9E-3 && ...
                abs(p2(1)) > 5E-6 && timeAboveDiff < 5E-3
            %}
            
            %plot(shapedPulse)
            %assignin('base','baseline',baseline);% - ...
            %mean(shapedPulse2(BLIndexStart-1000:BLIndexStart)) %+ ...
            %mean(shapedPulse(BLIndexEnd:BLIndexEnd+1000)))/2;% - ...
            %(mean(shapedPulse(index(1)-round(200E-9/T):index(1)-round(100E-9/T))) + ...
            %mean(shapedPulse(index(4)+round(100E-9/T):index(4)+round(200E-9/T))))/2;
            %{
            if amp >= 6
            assignin('base','p',p(1))
            evalin('base','ps = horzcat(ps,p);');
            end
            %}
            shapedPulse = shapedPulse/100;
            %plot(shapedPulse)
            %plot([pulseFit1(indexBegin:indexEnd) values(indexBegin:numel(x1))])
            %plot(shapedPulse)
            %plot(shapedPulse)
            %plot(shapedPulse(indexBegin:indexEnd));
            %indexBegin+round(1E-6/pulseIn.desc.Ts)
            %assignin('base','pulse',pulseIn.y);
            %evalin('base','pulses = pulses + pulse;');
            %+round(1.5E-6/pulseIn.desc.Ts)
            %{
            rng;
            x = rand;
            
            if x >= 0
                assignin('base','pulse',pulseIn.y(round(14E-6/T):round(17E-6/T)));
                evalin('base','pulses = horzcat(pulses,pulse);');
                
            end
            %}
            %{
        else
            shapedPulse = pulseIn.y;
            noGood = 1;
            amp = 0;
            %'BAD'
            %(difference(index(3)) - difference(index(2)))/(index(3) - index(2))
            %numel(index)
            
        end
        %}
        %plot(pulseDeconv.y)
        %shapedPulse = shapedPulse/100;
        %varianceUnDiff = mean(pulseDeconv.y(index(1)-round(4E-6/pulseIn.desc.Ts):index(1)-1000).^2)
        %varianceDiff = mean(shapedPulse(index(1)-round(4E-6/pulseIn.desc.Ts):index(1)-1000).^2)
        
        if plottingOn == 1
            %assignin('base','baseline',baseline);
            %evalin('base','baselines = horzcat(baselines,baseline);');
            %{
            assignin('base','p',pulseIn.y);
            evalin('base','ps =  horzcat(ps,p);');
            assignin('base','baseline',baseline);
            evalin('base','baseline = horzcat(baselines,baseline);');
            %}
            %assignin('base','shapedPulse',shapedPulse/100);
            %plot(shapedPulse/100,'-p','MarkerIndices',[index(1) index(2) index(3) index(4)], ...
            %    'MarkerFaceColor','red','MarkerSize',15);
            
            if amp >= 6.5
                assignin('base','pA',oldPulse(round(10E-6/T):round(35E-6/T)));
                evalin('base','pAs = horzcat(pAs,pA);');
            %{
            elseif amp >= 1.04736 && amp <= 1.08032
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
                %}
            end
            %}
            %{
                
            elseif amp > 4.3716 && amp <= 4.40917
                'AMP'
                amp
                assignin('base','p6',pulseIn.y(round(10E-6/T):round(27E-6/T)));
                evalin('base','p6s =  horzcat(p6s,p6);');
                
                
            elseif amp > 4.3716 && amp <= 4.40917
                'AMP'
                amp
                assignin('base','p7',pulseIn.y(round(10E-6/T):round(27E-6/T)));
                evalin('base','p7s =  horzcat(p7s,p7);');
                
            end
            %}
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