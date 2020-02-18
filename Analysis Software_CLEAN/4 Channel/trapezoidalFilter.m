function [trapezoidalPlateauAverage,vOut,noGood] = trapezoidalFilter(vIn,oldPulse,shapingTypeIn,TFRFIn,TFTopIn,tTiming)
format long
digits(15)

try
% --------------------------------------------------------------------------
% Trapezoidal filter for germanium pulses.
% If resolution gets worse, go back to a minimum value of peakDistance and
% only take peaks 1 and 2 in findPeaks.
% --------------------------------------------------------------------------

% Constants determining width of the flat top of the
% trapezoid (G) and the leading and falling edges of the
% trapezoid (L).
p = 1;%round(q/(L*pulseIn.desc.fs))
q = floor((1.25*p)/0.1);
Fs = (p/q)*vIn.desc.fs;
Ts = 1/Fs;
noGood = 0;

divid = 250E-9;

gRound = floor(TFTopIn/(divid));

L = round(TFRFIn/Ts);
G = round(gRound*divid/Ts);

[VMax,VMaxIndex] = max(vIn.y(50:numel(vIn.y)-50));
[VMin,VMinIndex] = min(vIn.y(50:numel(vIn.y)-50));

if abs(VMax) > abs(VMin) && tTiming == 0
    vIn.y = resample(vIn.y,p,q);
    oldPulse = resample(oldPulse,p,q);
    
else
    VMax = VMin;
    VMaxIndex = VMinIndex;
    
end

if abs(vIn.x(VMaxIndex)) >= 10E-6
    trapezoidalPlateauAverage = 0;
    vOut = [];
    noGood = 1;
    'Pulse not in coincidence'
    return;
    
end

offset = -vIn.info.OFFSET;
offset = offset + 250E-9;

N = numel(vIn.y);

V = abs(VMax);
VLower = 0.001;
VUpper = 0.5;

tau = (58000/(Ts/(10.^(-9))));

'T'
T
if T == 1
    
    % -------------------------------------------------------------------
    % Pole-zero cancelation
    % -------------------------------------------------------------------
    
    if shapingTypeIn == 50
       
        warnMsg = '';
            vInP.y = vIn.y;
        
        for n=2:numel(vIn.y)-1
            vInP.y(n) = vInP.y(n-1) + vIn.y(n) - vIn.y(n-1) + vIn.y(n-1)/tau;
            
        end
        
    else
        vInP.y = vIn.y;
        
    end
    
    % ---------------------------------------------------------------------
    % Recursive Trapezoid
    % ---------------------------------------------------------------------
    
    if shapingTypeIn == 1
        
        l1 = G + L;
        k1 = L;
        
        % Padding beginning and ending of vector with zeros to allow averaging "boxes" to
        % fit
        
        vInPNew.y = vertcat(zeros(2*L+G,1),vIn.y);
        endPadding = zeros(2*L+G,1);
        endPadding = endPadding + vInPNew.y(numel(vInPNew.y)-round((0.1*10.^(-6))/(Ts)));
        vInPNew.y = vertcat(vInPNew.y,endPadding);
        
        vOut = zeros(1,numel(vInPNew.y));
        vIn.y = vInPNew.y;
        
        d = zeros(1,numel(vInPNew.y));
        p = zeros(1,numel(vInPNew.y));
        r = zeros(1,numel(vInPNew.y));
        s = zeros(1,numel(vInPNew.y));
        for n=2:numel(vInPNew.y)-1
            
            if n-l1-k1 > 0
                d(n) = vIn.y(n) - vIn.y(n-k1) - vIn.y(n-l1) + vIn.y(n-k1-l1);
                p(n) = p(n-1) + d(n);
                r(n) = p(n) + tau*d(n);
                s(n) = s(n-1) + r(n);
                
            end
        end
        
        vOut = s*2/1E5;
        vOut = vOut(2*L+G:numel(vOut));
        vInPrint = vInP.y*(max(vOut)/max(vInP.y));
        vOutPrint = (vOut(1:numel(vInP.y)))';
        
        try
            baseline = (mean(vOut(newsIndex(1)-500:newsIndex(1)-100)) + mean(vOut(newsIndex(4)+100:newsIndex(4)+500)))/2;
            trapezoidalPlateauAverage = mean(vOut(newsIndex(3):newsIndex(3)));
            vOut = horzcat(vInPrint,vOutPrint);
        catch
            trapezoidalPlateauAverage = 0;
            vOut = horzcat(vInPrint,vOutPrint);
        end
        
        % --------------------------------------------------------------------
        % Integration
        % --------------------------------------------------------------------
        
    elseif shapingTypeIn == 0
        
        netInput = oldPulse(round(10E-6/Ts):round(35E-6/Ts));
        baseline = mean(netInput(200:450));
        netInput = netInput - baseline;
        
        netInput = netInput(400:1200);
        
        plot(netInput)
        net = evalin('base','net3');
        
        out = net(netInput);
        assignin('base','out',out);
        assignin('base','netInput',netInput);
        
        out = smooth(out*100,50,'moving');
        
        trapezoidalPlateauAverage = max(out(1:500));
        vOut = out;
        noGood = 0;
        
        if trapezoidalPlateauAverage <= 1.1
            trapezoidalPlateauAverage = 0;
            
        end
        
        % --------------------------------------------------------------------
        % Manual convolution
        % --------------------------------------------------------------------
    elseif shapingTypeIn == 2
        
        try
            measurePulse = smooth(vIn.y,101,'moving') - ...
                mean(vIn.y(round(offset/Ts-2E-6/Ts):round(offset/Ts-1E-6/Ts)));
            
            [m,ending] = max(measurePulse);
            [~,begin] = min(abs(measurePulse(round(offset/Ts-1E-6/Ts):round(offset/Ts))-0.1*m));
            
            begin = begin + round((offset/Ts-1E-6/Ts));
            
            riseTime = (vIn.x(round(ending*q/p)) - vIn.x(round(begin*q/p)));
            
            if riseTime > 1100E-9 || riseTime < 110E-9
                noGood = 1;
                trapezoidalPlateauAverage = 0;
                vOut = vIn.y;
                'rise time'
                return;
                
            end
            
        catch
            noGood = 1;
            trapezoidalPlateauAverage = 0;
            vOut = vIn.y;
            'rise time error'
            return;
        end
        
        pulse = vIn.y(round(11E-6/Ts):round(27E-6/Ts))' + 1;
        
        r = evalin('base','r');
        undone = log(pulse);
        r = r*max(undone);
        
        undone(550:numel(undone)) = undone(550:numel(undone)) - r';
        
        maxIndex = round(offset/Ts-round(9.5E-6/Ts));
        endLinear = round(10E-6/Ts);
        undo2 = undone(maxIndex:maxIndex+endLinear-1);
        x = linspace(maxIndex,maxIndex+endLinear,endLinear);
        %assignin('base','x',x)
        %undone = undo2;
        
        p = polyfit(x,real(undo2),1);
        assignin('base','undo2',undo2)
        assignin('base','x',x)
        
        pVal = polyval(p,x);%feval(fitresult,x);%polyval(p,x);
        
        trapezoidalPlateauAverage = pVal(1)*100 - mean(undone(200:350))*100
        vOut = undone;
        
        % ---------------------------------------------------------------------
        % Moving Average Semi-Gaussian
        % ---------------------------------------------------------------------
        
    elseif shapingTypeIn == 3
        
        vOut = abs(smooth(diff(vIn.y),L,'moving'));
        
        for l=1:10
            vOut = smooth(vOut(100:numel(vOut)-100),L,'moving');
        end
        
        vOut = vOut(100:numel(vOut)-3000)*1000*(L);
        [trapezoidalPlateauAverage,index] = max(vOut);
        trapezoidalPlateauAverage
        
        if  isempty(index) == 1 || isempty(trapezoidalPlateauAverage) == 1 || ...
                index == 1
            trapezoidalPlateauAverage = 0;
            
        end
        assignin('base','vIn',abs(vIn.y));
        assignin('base','vOut',vOut/(1000*(L)));
        
        % --------------------------------------------------------------------
        % Moving window deconvolution.
        % --------------------------------------------------------------------
        
    elseif shapingTypeIn == 4
        
        try
            measurePulse = smooth(vIn.y,101,'moving') - ...
                mean(vIn.y(round(offset/Ts-2E-6/Ts):round(offset/Ts-1E-6/Ts)));
            
            [m,ending] = max(measurePulse);
            [~,begin] = min(abs(measurePulse(round(offset/Ts-1E-6/Ts):round(offset/Ts))-0.1*m));
            
            begin = begin + round((offset/Ts-1E-6/Ts));
            
            riseTime = (vIn.x(round(ending*q/p)) - vIn.x(round(begin*q/p)))
            
            if riseTime > 1100E-9 || riseTime < 110E-9
                noGood = 1;
                trapezoidalPlateauAverage = 0;
                vOut = vIn.y;
                'rise time'
                return;
                
            end
            
        catch
            noGood = 1;
            trapezoidalPlateauAverage = 0;
            vOut = vIn.y;
            'rise time error'
            return;
        end
        
        ap = 0;
        newvIn = vIn.y(round(12E-6/Ts):round(24E-6/Ts));
        for n=1:numel(newvIn)-1
            assignin('base','pulse',newvIn);
            assignin('base','n',n);
            ap = ap + evalin('base','h(n)*pulse(n);');
            
        end
        
        assignin('base','ap',ap);
        newestvIn = real(evalin('base','ap.*s;'));
        newestvIn2 = padarray(newestvIn,round(12E-6/Ts)-1,0,'pre');
        
        [ml,m] = max(newestvIn);
        
        trapezoidalPlateauAverage = ml*100;
        vOut = [newvIn newestvIn'];
        
    elseif shapingTypeIn == 5
        
         try
            measurePulse = smooth(vIn.y,101,'moving') - ...
                mean(vIn.y(round(offset/Ts-2E-6/Ts):round(offset/Ts-1E-6/Ts)));
            
            [m,ending] = max(measurePulse);
            [~,begin] = min(abs(measurePulse(round(offset/Ts-1E-6/Ts):round(offset/Ts))-0.1*m));
            
            begin = begin + round((offset/Ts-1E-6/Ts));
            
            riseTime = (vIn.x(round(ending*q/p)) - vIn.x(round(begin*q/p)));
            
            if riseTime > 1100E-9 || riseTime < 110E-9
                noGood = 1;
                trapezoidalPlateauAverage = 0;
                vOut = vIn.y;
                'rise time'
                return;
                
            end
            
        catch
            noGood = 1;
            trapezoidalPlateauAverage = 0;
            vOut = vIn.y;
            'rise time error'
            return;
        end
        
        try
            baseline = mean(vIn.y(round((offset-999E-9)*Fs):round((offset-500E-9)*Fs)));
            vIn.y = vIn.y;
            
            [vIn.y,err,a] = MLSmoothing(vIn.y,Ts);
            err
            
            if err < 2.5E-3
                trapezoidalPlateauAverage = a*100;
                vOut = vIn.y*100;
                assignin('base','vOut',vOut/100);
                
            else
                vOut = [];
                
                trapezoidalPlateauAverage = 0;
                
                noGood = 1;
                
            end
            
        catch
            vOut = [];
            
            trapezoidalPlateauAverage = 0;
            
            noGood = 1;
            
            'MLSmoothing error'
            
        end
        
        % -------------------------------------------------------------------
        % Fitting function
        % -------------------------------------------------------------------
        
    elseif shapingTypeIn == 7
        
        xAxis = linspace(1,numel(vIn.y),numel(vIn.y))';
        
        v = vIn.y(round(10E-6/Ts):round(24E-6/Ts));
        
        [fitResult, gof] = fitting3(xAxis(round(10E-6/Ts):round(24E-6/Ts)),v);
        gof.adjrsquare
        
        double a;
        a = double(vpa(coeffvalues(fitResult)));
        assignin('base','a',a);
        y = feval(fitResult,xAxis);
        amp = (2*a(1))*100;
        %vIns = evalin('base','vIn;');
        
        if abs(gof.adjrsquare) >= 0.998
            
            trapezoidalPlateauAverage = amp;
            vOut = [v*10.^3,feval(fitResult,xAxis(round(10E-6/Ts):round(24E-6/Ts)))*10.^3];
            
        else
            
            vOut = [];
            trapezoidalPlateauAverage = 0;
            
        end
        
        % Parabola fitting for shaped pulse.
        
    elseif shapingTypeIn == 9
        
        [~,pulseMaxIndex] = max(vIn.y);
        
        pulseRangeMin = pulseMaxIndex - L;
        pulseRangeMax = pulseMaxIndex + L;
        
        if pulseRangeMax < numel(vIn.y) && pulseRangeMin > 0
            
            lastwarn('')
            
            [p,S,mu] = polyfit(vIn.x(pulseRangeMin:pulseRangeMax), ...
                vIn.y(pulseRangeMin:pulseRangeMax),2);
            assignin('base','p',p);
            [warnmsg,~] = lastwarn;
            
            [vOut,delta] = polyval(p,vIn.x(pulseRangeMin:pulseRangeMax),S,mu);
            
            vOut = [vIn.y(pulseRangeMin:pulseRangeMax),vOut]*100;
            
            RSS = sum(delta.^2);
            
            baseline = mean(vIn.y(100:100+round(10E-6/Ts)));
            
            if isempty(warnmsg) && RSS <= 2.9*10.^(-4)
                values = p;
                
                trapezoidalPlateauAverage = ((values(3) - values(2)^2/(4*values(1))))*100;
                
            else
                trapezoidalPlateauAverage = 0;
                
            end
            
        else
            trapezoidalPlateauAverage = 0;
            vOut = [];
            
        end
        
    elseif shapingTypeIn == 10
        try
            measurePulse = smooth(vIn.y,101,'moving') - ...
                mean(vIn.y(round(offset/Ts-2E-6/Ts):round(offset/Ts-1E-6/Ts)));
            
            [m,ending] = max(measurePulse);
            [~,begin] = min(abs(measurePulse(round(offset/Ts-1E-6/Ts):round(offset/Ts))-0.1*m));
            
            begin = begin + round((offset/Ts-1E-6/Ts));
            
            riseTime = (vIn.x(round(ending*q/p)) - vIn.x(round(begin*q/p)));
            
            if riseTime > 1100E-9 || riseTime < 110E-9
                noGood = 1;
                trapezoidalPlateauAverage = 0;
                vOut = vIn.y;
                'rise time'
                return;
                
            end
            
        catch
            noGood = 1;
            trapezoidalPlateauAverage = 0;
            vOut = vIn.y;
            'rise time error'
            return;
        end
        
        baseline = mean(vIn.y(round((offset-999E-9)*Fs):round((offset-500E-9)*Fs)));
        
        v = vIn.y(round(10E-6/Ts)+240:round(30E-6/Ts)-1084) - baseline;
        size(v)
        
        assignin('base','vIns',vIn.y(round(10E-6/Ts):round(27E-6/Ts)-771));
        assignin('base','vIns2',v);
        assignin('base','vIns3',vIn.y(round(10E-6/Ts)+399:round(27E-6/Ts)-771));
        aL = evalin('base','net(vIns2);');
        assignin('base','aL',aL);
        
        plot([v aL])
        err = immse(aL,v)
        
        aF = aL; 
        
        if err <= 1E-6
        trapezoidalPlateauAverage = max(aF)*100
        
        amp = trapezoidalPlateauAverage;
        vOut = aL;
        
        else
            trapezoidalPlateauAverage = 0;
            vOut = v;
            err
            
        end
        
        if false
        if amp >= 0.9297 && amp <= 0.9519
                'AMP'
                amp
                %figure
                %plot(pulseIn.y);
                %title('low side');
                assignin('base','p1',vIn.y(round(10E-6/Ts):round(24E-6/Ts)));
                evalin('base','p1s = horzcat(p1s,p1);');
                
            elseif amp > 3.185 && amp <= 3.2078
                'AMP'
                amp
                assignin('base','p2',vIn.y(round(10E-6/Ts):round(24E-6/Ts)));
                evalin('base','p2s =  horzcat(p2s,p2);');
                
            elseif amp >= 3.5137 && amp <= 3.5461
                'AMP'
                amp
                %assignin('base','pulse',pulseIn.y);
                %evalin('base','pulseHigh = pulseHigh + pulse;');
                %evalin('base','m=m+1');
                assignin('base','p3',vIn.y(round(10E-6/Ts):round(24E-6/Ts)));
                evalin('base','p3s =  horzcat(p3s,p3);');
                
            elseif amp > 4.1306 && amp <= 4.1648
                'AMP'
                amp
                assignin('base','p4',vIn.y(round(10E-6/Ts):round(24E-6/Ts)));
                evalin('base','p4s =  horzcat(p4s,p4);');
                
            elseif amp > 4.428 && amp <= 4.4536
                'AMP'
                amp
                assignin('base','p5',vIn.y(round(10E-6/Ts):round(24E-6/Ts)));
                evalin('base','p5s =  horzcat(p5s,p5);');
                
            end
        end
        
    end
    
else
    trapezoidalPlateauAverage = 0;
    vOut = [];
    
end

catch
    trapezoidalPlateauAverage = 0;
    vOut = [];
    noGood = 1;
    'trapezoidal filter error'
    
end

end