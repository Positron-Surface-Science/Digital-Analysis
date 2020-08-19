function [trapezoidalPlateauAverage,vOut,noGood] = trapezoidalFilter(vIn,oldPulse,shapingTypeIn,TFRFIn,TFTopIn,tTiming)
format long
digits(15)

%try
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

L = round(TFRFIn/Ts); %3*10^(-6)
G = round(gRound*divid/Ts);

%vIn.y = smooth(vIn.y,Ts/vIn.desc.fs,'moving');

[VMax,VMaxIndex] = max(vIn.y(50:numel(vIn.y)-50));
[VMin,VMinIndex] = min(vIn.y(50:numel(vIn.y)-50));

if ~isempty(VMax) && ~isempty(VMin) && abs(VMax) > abs(VMin) && tTiming == 0
    vIn.y = resample(vIn.y,p,q);
    oldPulse = resample(oldPulse,p,q);
    %vIn.y = resample(vIn.y,p,q);
    
elseif ~isempty(VMax) && ~isempty(VMin)
    VMax = VMin;
    VMaxIndex = VMinIndex;
    
else
    trapezoidalPlateauAverage = 0;
    vOut = [];
    noGood = true;
    return;
    
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

tau = (58000/(Ts/(10.^(-9))));%29250; % Decay constant of HPGe preamplifier. 58000
%tau
% Requiring a minimum and maximum voltage to reduce noise and
% over-pulses
%try
    T = 1;%V > VLower && V < VUpper;
%catch
    %trapezoidalPlateauAverage = 0;
    %vOut = [];
    %noGood = 1;
    
%end
'T'
T
if T == 1
    
    % -------------------------------------------------------------------
    % Pole-zero cancelation
    % -------------------------------------------------------------------
    
    if shapingTypeIn == 50
        %{
        xAxis = linspace(1,numel(vIn.y),numel(vIn.y))';
        %vIn.y = smooth(vIn.y,5001,'lowess');
        
        [fitResult, gof] = fitting(xAxis,vIn.y);
        gof.adjrsquare
        
        double a;
        a = coeffvalues(fitResult);
        assignin('base','a',a);
        %}
        %if gof.adjrsquare >= 0.99
        
        %   vFit = [vIn.y*10.^3,feval(fitResult,xAxis)*10.^3];
        
        %end
        
        warnMsg = '';
            vInP.y = vIn.y;
        
        for n=2:numel(vIn.y)-1
            %sumGz = sumGz + vIn.y(n) - vIn.y(n-1) + (vIn.y(n-1)/tau);
            vInP.y(n) = vInP.y(n-1) + vIn.y(n) - vIn.y(n-1) + vIn.y(n-1)/tau;
            
        end
        %{
        figure
        plot(vOut)
        figure
        plot(vOutDiff2);
        %}
        %{
        try
            trapezoidalMax = max(abs(vOutDiff2(1:numel(vOutDiff2)-1000)));
            
            halfTrapezoidalMax = trapezoidalMax/2;
            
            [~,newsIndex] = findpeaks(abs(vOutDiff2(1:numel(vOutDiff2)-1000)), ...
                'MinPeakHeight',halfTrapezoidalMax,'MinPeakDistance', ...
                25E-9/vIn.desc.Ts);
            %}
            % Removing pulses on a distorted baseline
            %{
            beforePulseBaseline = mean(vIn.y(round(newsIndex(1)- ...
                (5E-6/vIn.desc.Ts):newsIndex(1))));
            
            if beforePulseBaseline >= 0.01
                vOut = vIn.y;
                trapezoidalPlateauAverage = 0;
                'before'
                return;
                
            end
            
            if numel(newsIndex) == 2
                newsIndexBegin = newsIndex(2);%+(500*10.^(-9)/vIn.desc.Ts));
                newsIndexEnd = numel(vIn.y) - 100;
                
            elseif numel(newsIndex) == 4
                newsIndexBegin = newsIndex(2);%+(500*10.^(-9)/vIn.desc.Ts));
                newsIndexEnd = newsIndex(3);
                'pileup'
                
            else
                vOut = vIn.y;
                trapezoidalPlateauAverage = 0;
                'peaks'
                return;
                
            end
            %}
            %{
            newsIndexBegin = newsIndex(2);
            newsIndexEnd = numel(vInP.y) - 100;
            %}
            %{
            % ----------------------------------------------------------
            % Probably bad
            % ----------------------------------------------------------
            
            baselineRange = round(newsIndex(1)-(1E-6/vIn.desc.Ts));
            vOrig = vIn.y;
            lastwarn('')
            baselineFit = polyfit(vIn.x(100:baselineRange), ...
                vIn.y(100:baselineRange),1);
            warnMsg = lastwarn;
            
            if baselineFit(1) < -10
                baseline = polyval(baselineFit,vIn.x);
                %vInPUncorrected = vIn.y;
                vIn.y = vIn.y - baseline;
                
            end
            %assignin('base','baseline',baseline)
            
            % ----------------------------------------------------------
            % Probably bad end
            % ----------------------------------------------------------
            %}
            %plot(vIn.x,vIn.y,vIn.x,vOrig,vIn.x,baseline)
            %tau = a(4);
            
        
            %baseline = mean(vOut(numel(vOut)-5000): ...
            %    numel(vOut)-1000);
            
            %vOut = (vOut - baseline);
            %}
            %{
            lastwarn('')
            pFit = polyfit(vIn.x(newsIndexBegin: ...
                newsIndexEnd),vInP.y(newsIndexBegin:newsIndexEnd),1);
            warnMsg2 = lastwarn;
            %}
            %{
            %'before correction'
            %pFit(1)
            
            %p = polyval(pFit,vIn.x);
            %assignin('base','p',p);
            %figure
            %plot(vOutCorrection(:,1));
            %}
            % Applying a correction to the time constant of the pulse to
            % improve the exponential deconvolution.
            %{
            vInPOriginal = vInP.y;
            tauCounter = 1;
            number = 500;
            
            while (pFit(1) > 1 || pFit(1) < -1) && tauCounter <= 20
                vInP.y = vInPOriginal;
                if pFit(1) > 1
                    tau = tau + (number/(vIn.desc.Ts/(10.^(-9))));
                elseif pFit(1) < -1
                    tau = tau - (number/(vIn.desc.Ts/(10.^(-9))));
                end
                
                for n=2:numel(vIn.y)-1
                    %sumGz = sumGz + vIn.y(n) - vIn.y(n-1) + (vIn.y(n-1)/tau);
                    vInP.y(n) = vInP.y(n-1) + vIn.y(n) - vIn.y(n-1) +
                    vIn.y(n-1)/tau;
                    
                end
                
                %figure
                %plot(vOutCorrection(:,tauCounter));
                lastwarn('')
                pFit = polyfit(vIn.x(newsIndexBegin: ...
                    newsIndexEnd),vInP.y(newsIndexBegin:newsIndexEnd),1);
                warnMsg2 = lastwarn;
                
                %p = polyval(pFit,vIn.x);
                tauCounter = tauCounter + 1;
                
            end
            %{
            %plot(vIn.x,vInPOriginal,vIn.x,vInP.y,'-p','MarkerIndices',[newsIndexBegin newsIndexEnd],'MarkerFaceColor','red','MarkerSize',15)
            %{
            pFit(1)
            tauCounter
            plot(vInP.y)
            %}
        %plot(vIn.x,vInPUncorrected,vIn.x,vInP.y,vIn.x,baseline);
       %}
        catch
            vOut = vIn.y;
            trapezoidalPlateauAverage = 0;
            'error'
            return;
            
        end
        
        %{
        if isempty(warnMsg) == 0 || isempty(warnMsg2) == 0
            vOut = vInP.y;
            trapezoidalPlateauAverage = 0;
            'warning'
            return;
            
        end
        %}
        if tauCounter >= 20
            vOut = vInP.y;
            trapezoidalPlateauAverage = 0;
            'distorted'
            return;
            
        end
        
        %{
        'after correction'
        pFit(1)
        tauCounter
        figure
        
        assignin('base','vInx',vIn.x);
        assignin('base','p',p);
        assignin('base','vInP',vInP.y);
        %}
        %figure
        %plot(vIn.x,p,vIn.x,vInP.y);
        %figure
        %plot(vIn.x,vInPOriginal)
        
        
        %figure;
        %plot(vInP.y);
        
        %vInP.y = smooth(vInP.y,5001,'sgolay',9);
        %vInP.y = smooth(vInP.y,201,'sgolay',2);
        %vInP.y = smooth(vInP.y,51,'sgolay',1);
        
        %figure;
        %plot(vInP.y);
       %}
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
        
        %vInP.y = vertcat(vInP.y,zeros(2*L+G,1));
        
        vOut = zeros(1,numel(vInPNew.y));
        vIn.y = vInPNew.y;
        %figure;
        %plot(vInP.y);
        d = zeros(1,numel(vInPNew.y));
        p = zeros(1,numel(vInPNew.y));
        r = zeros(1,numel(vInPNew.y));
        s = zeros(1,numel(vInPNew.y));
        for n=2:numel(vInPNew.y)-1
            %numel(vIn.y)
            %n-1
            %n-k1
            %n-l1
            %n-l1-k1
            if n-l1-k1 > 0
                %vOut(n) = vOut(n-1) + vInP.y(n) - vInP.y(n-k1) - vInP.y(n-l1-k1);
                %vOut(n) = vOut(n-1) + (vInPNew.y(n) - vInPNew.y(n-k1)) - (vInPNew.y(n-l1) - vInPNew.y(n-l1-k1));
                d(n) = vIn.y(n) - vIn.y(n-k1) - vIn.y(n-l1) + vIn.y(n-k1-l1);
                p(n) = p(n-1) + d(n);
                r(n) = p(n) + tau*d(n);
                s(n) = s(n-1) + r(n);
                
            end
        end
        
        %figure
        %plot(vOut)
        
        %vOut = smooth(vOut,1051,'sgolay',2);
        %figure
        %plot(vOut)
        
        vOut = s*2/1E5;
        
        % Correcting the baseline after shaping
        
        %baseline = mean(vOut(1:500));
        
        %vOut = vOut - baseline;
        
        %trapezoidalPlateauAverage = max(vOut);
        
        vOut = vOut(2*L+G:numel(vOut));
        vInPrint = vInP.y*(max(vOut)/max(vInP.y));
        vOutPrint = (vOut(1:numel(vInP.y)))';
        %assignin('base','vInPrint',vInPrint);
        %assignin('base','vOutPrint',vOutPrint);
        
        %{
        
        vOutDiff1 = diff(vOut);
        vOutDiff1Smooth = smooth(vOutDiff1,501,'lowess');
        vOutDiff2 = diff(vOutDiff1Smooth);
        
        trapezoidalMax = max(abs(vOutDiff2(1:numel(vOutDiff2)-1000)));
        
        halfTrapezoidalMax = trapezoidalMax/2;
        
        [~,newsIndex] = findpeaks(abs(vOutDiff2(1:numel(vOutDiff2)-1000)), ...
            'MinPeakHeight',halfTrapezoidalMax);
        %}
        try
            baseline = (mean(vOut(newsIndex(1)-500:newsIndex(1)-100)) + mean(vOut(newsIndex(4)+100:newsIndex(4)+500)))/2;
            trapezoidalPlateauAverage = mean(vOut(newsIndex(3):newsIndex(3)));
            vOut = horzcat(vInPrint,vOutPrint);
        catch
            trapezoidalPlateauAverage = 0;
            vOut = horzcat(vInPrint,vOutPrint);
        end
        
        %plot(vOut)
    %{
    %figure;
    %plot(abs(vOutDiff2(1:numel(vOutDiff2)-1000)));
    
    try
        trapezoidalMax = max(abs(vOutDiff2(1:numel(vOutDiff2)-1000)));
        
        halfTrapezoidalMax = trapezoidalMax/2;
        
        [~,newsIndex] = findpeaks(abs(vOutDiff2(1:numel(vOutDiff2)-1000)), ...
            'MinPeakDistance',0.9*G,'MinPeakHeight',halfTrapezoidalMax);
        
        baseline = mean(vOut(1:1000));
        
        vOut = vOut - baseline;
        
        trapezoidalPlateauAverage = mean(vOut(newsIndex(2):newsIndex(3)));
        
    catch
        trapezoidalPlateauAverage = 0;
        
    end
    
        %}
        %{
    numel(vIn.x)
    numel(vOut)
    figure;
    plot(vIn.x,vOut,'-p','MarkerIndices',[newsIndex(2) newsIndex(3)], ...
    'MarkerFaceColor','red','MarkerSize',15);
    
    %assignin('base','news',news);
    assignin('base','newsIndex',newsIndex);
    assignin('base','vOut',vOut);
        %}
        
        %{
if isempty(squarePulse)

    squarePulse = zeros(numel(vIn.y),1);

    for i=1:numel(vIn.y)
        if i <= L
            squarePulse(i) = -1;
        elseif i <= L + G
            squarePulse(i) = 0;
        elseif i <= 2*L + G
            squarePulse(i) = 1;
        else
            squarePulse(i) = 0;
        end
    end
end

vOut = conv(vIn.y(50:numel(vIn.y)-50),squarePulse);
        %}
        
        % --------------------------------------------------------------------
        % Integration
        % --------------------------------------------------------------------
        
    elseif shapingTypeIn == 0
        
        %{
        vInDifferentiated = vIn.y;
        vInSmoothed = smooth(vInDifferentiated,1001,'lowess');
        for j=1:15
            vInSmoothed = smooth(vInSmoothed,1001,'lowess');
        end
        vInSmoothed = vInSmoothed - mean(vInSmoothed(100:5000));
        vOut = max(vInSmoothed);%abs(trapz(vIn.x,vIn.y)/vIn.desc.Ts);
        %figure
        %plot(vInSmoothed)
        trapezoidalPlateauAverage = vOut;
        %}
        
        %[vIn.y,~,~] = MLSmoothing(vIn.y,Ts);
        %{
        netInput = oldPulse(round(10E-6/Ts):round(35E-6/Ts));
        baseline = mean(netInput(200:450));
        netInput = netInput - baseline;
        %netInput = netInput/max(netInput);
        netInput = netInput(400:1200);
        %netInput = netInput(1:50);
        
        %netInput = netInput(15:40);
        plot(netInput)
        net = evalin('base','net3');
        
        out = net(netInput);
        assignin('base','out',out);
        assignin('base','netInput',netInput);
        %}
        %truth = [1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1 1 1 0 1 0];
        %truth = evalin('base','truth');
        %result = truth*out;
        %{
        out = smooth(out*100,50,'moving');
        
        trapezoidalPlateauAverage = max(out(1:500));
        vOut = out;
        noGood = 0;
        
        if trapezoidalPlateauAverage <= 1.1
            trapezoidalPlateauAverage = 0;
            
        end
        %}
        %if sum(out(1:30)) == 1 %out(8) == 1 || out(13) == 1 || out(17) == 1 || out(22) == 1
        %plot(vIn.y)
        %}
        [trapezoidalPlateauAverage,vOut,noGood] = MWD(vIn,oldPulse,TFRFIn,TFTopIn);
        
        %{
        else
            
            trapezoidalPlateauAverage = 0;
            vOut = vIn.y;
            noGood = 1;
            
        end
        %}
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
        
        %[fitresult,gof] = createFit(x,undo2);
        
        pVal = polyval(p,x);%feval(fitresult,x);%polyval(p,x);
        
        trapezoidalPlateauAverage = pVal(1)*100 - mean(undone(200:350))*100
        vOut = undone;
        
        %{
        vInP.y = vertcat(zeros(2*L+G,1),vInP.y);
        
        v = vInP.y;
        
        Vav1 = 0;
        Vav2 = 0;
        
        % Algorithm for building the trapezoid (calculating a
        % moving average).
        
        % Do not attempt to calculate average when leading box
        % reaches the end of the pulse.
        for n=1:numel(vIn.y)-L
            
            Vinav1 = 0;
            Vinav2 = 0;
            
            % Begin averaging of the leading box before trailing
            % box enters the pulse range (trailing box average is
            % 0).
            if n <= G
                for j=1:L
                    Vinav2 = Vinav2 + v(n+j);
                    Vinav1 = 0;
                end
                
                % Continue averaging the leading box while beginning to
                % average the fractional part of the trailing box
                % (average of thepart of the trailing box outside of
                % the range of the pulse is 0).
            elseif n >= G && n - G <= L
                for j=1:L
                    Vinav2 = Vinav2 + v(n+j);
                end
                for j=1:n-G
                    Vinav1 = Vinav1 + v(j);
                end
                
                % Continue average with both averaging boxes inside of
                % the range of the pulse.
            elseif n >= G && n - G >= L
                for j=1:L
                    Vinav2 = Vinav2 + v(n+j);
                    Vinav1 = Vinav1 + v(n-G-L+j);
                end
                
            end
            
            % Averaging over L values.
            Vav1(n) = Vinav1/L;
            Vav2(n) = Vinav2/L;
            
            % Subtracting averages to create trapezoid.
            vOut(n) = Vav2(n) - Vav1(n);
            
        end
        
        vOut = vOut - mean(vOut(1:1000));
        
        vOut = vOut*10.^4;
        trapezoidalPlateauAverage = max(vOut);
        %}
        % ---------------------------------------------------------------------
        % Moving Average Semi-Gaussian
        % ---------------------------------------------------------------------
        
    elseif shapingTypeIn == 3
        %vIn.y = vIn.y - mean(vIn.y(100:1000));
        %Fs = vIn.desc.fs; % Sampling frequency.
        %N = numel(vInP.y); % Number of samples in waveform.
        %dN = N/Fs; % Time per channel.
        %timeDelay = L/dN % Time delay for differentiation.
        %vInP.y = vIn.y;
        %
        %{
        vInPNew.y = vertcat(zeros(2*L,1),vInP.y(1:numel(vInP.y)-10));
        endPadding = zeros(2*L+G,1);
        endPadding = endPadding + vInPNew.y(numel(vInPNew.y)-round((0.1*10.^(-6))/(vIn.desc.Ts)));
        vInPNew.y = vertcat(vInPNew.y,endPadding);
        diffGauss = vInPNew.y;
        %}
        %{
        vInP.y = vertcat(zeros(2*L,1),vInP.y);
        endPadding = zeros(2*L+G,1);
        endPadding = endPadding + vInP.y(numel(vInP.y-100));
        vInP.y = vertcat(vInP.y,endPadding);
        diffGauss = vInP.y;
        %}
        %{
        for i=L+1:numel(vInPNew.y)-10 %i=2*L+1:numel(vInP.y)-(2*L+G)
            diffGauss(i) = vInPNew.y(i) - vInPNew.y(i-L);%diff(vInP.y(100:numel(vInP.y)-100));
            
        end
        diffGauss = diffGauss(1:numel(diffGauss)-10);
        %}
        %assignin('base','diffGauss',diffGauss);
        
        %----------------------------------------------------------------------
        % Multi-Site Event Cancellation
        %----------------------------------------------------------------------
        %{
        baseline = mean(diffGauss(100:1000));
        [~,deltaIndex] = max(diffGauss(100:numel(diffGauss)-100));
        deltaIndex = deltaIndex + 100;
        sumLeft = 0;
        sumRight = 0;
        foundBaselineLeft = false;
        foundBaselineRight = false;
        i = 1;
        %}
        %figure;
        %plot(diffGauss);
        %{
        while i < (300/vIn.desc.Ts) && deltaIndex - i > 0 && deltaIndex + i < numel(diffGauss)
            
            if diffGauss(deltaIndex-i) > baseline && foundBaselineLeft == false
                sumLeft = sumLeft + diffGauss(deltaIndex-i);
                diffGauss(deltaIndex-i) = 0;
                
            elseif diffGauss(deltaIndex-i) <= baseline
                foundBaselineLeft = true;
            end
            
            if diffGauss(deltaIndex+i) > baseline && foundBaselineRight == false
                sumRight = sumRight + diffGauss(deltaIndex+i);
                diffGauss(deltaIndex+i) = 0;
                
            elseif diffGauss(deltaIndex+i) <= baseline
                foundBaselineRight = true;
            end
            
            i = i + 1;
            
        end
        
        sumPeak = sumLeft + sumRight;
        diffGauss(deltaIndex) = diffGauss(deltaIndex) + sumPeak;
        diffGaussX = linspace(1,numel(diffGauss),numel(diffGauss))';
        %}
        %figure;
        %plot(diffGauss);
        
        %{
    figure;
    plot(vInP.y(100:numel(vInP.y)-100));
    figure;
    plot(diffGauss);
        %}
        %{
        for i=1:15
            
            diffGauss = smooth(diffGauss,L,'moving');
            
            %{
        for j=2*L:numel(diffGauss)
            for k=j-L:j
                diffGauss(j) = diffGauss(k);
                
            end
            
            diffGauss = diffGauss/L;
            
        end
            %}
            %{
        for j=L/2+1:(numel(diffGauss)-L-1)
            diffGauss(i,j) = 0;
            
            for s=-L/2:L/2
                diffGauss(i,j) = diffGauss(i,j) + diffGauss(i-1,j+s);

            end
            
            diffGauss(i,j) = diffGauss(i,j)/L;
            
        end
            %}
        end
        
        %figure;
        %plot(diffGauss);
        
        baseline = mean(diffGauss(100:5000));
        diffGauss = diffGauss - baseline;
        %{
        [~,pulseMaxIndex] = max(diffGauss);
        
        pulseRangeMin = pulseMaxIndex - round((2.5*10.^(-6))/(vIn.desc.Ts));
        pulseRangeMax = pulseMaxIndex + round((2.5*10.^(-6))/(vIn.desc.Ts));
        
        %numel(diffGaussX(pulseRangeMin:pulseRangeMax))
        %numel(diffGauss(pulseRangeMin:pulseRangeMax))
        
        if pulseRangeMax < numel(diffGauss) && pulseRangeMin > 0
            
            lastwarn('')
            
            [p,~,~] = polyfit(diffGauss(pulseRangeMin:pulseRangeMax), ...
                diffGauss(pulseRangeMin:pulseRangeMax),2);
            
            [warnmsg,~] = lastwarn;
            
            if isempty(warnmsg)
                values = p;%coeffvalues(parabolaFit);
                
                trapezoidalPlateauAverage = (values(3) - values(2)^2/(4*values(1)))*100;
                
            else
                trapezoidalPlateauAverage = 0;
                
            end
            
        else
            trapezoidalPlateauAverage = 0;
            
        end
        %}
        
        vOut = (diffGauss(100:numel(diffGauss)-100)*10.^4)';
        
        trapezoidalPlateauAverage = max(vOut(1:numel(vOut)-1000));
        
        %figure;
        %plot(vOut);
        %}
        
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
        
        trapezoidalPlateauAverage = ml*100;% - mean(newvIn(100:275))*100;
        vOut = [newvIn newestvIn'];
        
        %plot(vOut)
        
        
        %{
        %diffGauss = diff(vInP.y(100:numel(vInP.y)-100));
        %vInP.y = vertcat(zeros(2*L,1),vInP.y);
        %{
        for i=L+1:numel(vInP.y)-10
            diffGauss(i) = vInP.y(i) - vInP.y(i-L);%diff(vInP.y(100:numel(vInP.y)-100));
            
        end
        %}%{
        movingSum = movsum(vIn.y,L)';
        for i=1:numel(vIn.y) %i=2*L+1:numel(vInP.y)-(2*L+G)
            if i > L + 1
                difference(i) = vIn.y(i) - vIn.y(i-L);%diff(vInP.y(100:numel(vInP.y)-100));
                vOut(i) = difference(i) + movingSum(i) - movingSum(i)/tau;
            else
                difference(i) = 0;
                vOut(i) = 0;
            end
            
        end
        %}
        
        %vInP.y = vIn.y;
        difference = vInP.y;
        vOut = vInP.y;
        
        for n=2:numel(vIn.y)-1
            %sumGz = sumGz + vIn.y(n) - vIn.y(n-1) + (vIn.y(n-1)/tau);
            %vInP.y(n) = vInP.y(n-1) + vIn.y(n) - vIn.y(n-1) + vIn.y(n-1)/tau;
            if n > G
                difference(n) = vInP.y(n) - vInP.y(n-G);
            else
                difference(n) = 0;
            end
            
            vOut(n) = difference(n)*100;
        end
        
        %assignin('base','vOut',vOut);
        %vOut = vInP.y;
        vOutDiff1 = diff(vOut);
        vOutDiff1Smooth = smooth(vOutDiff1,501,'lowess');
        vOutDiff2 = diff(vOutDiff1Smooth);
        
        %figure
        %plot(vOut)
        %figure
        
        %plot(vIn.x(1:numel(vIn.x)-1),abs(vOutDiff1Smooth));
        
        
        try
            trapezoidalMax = max(abs(vOutDiff2(1:numel(vOutDiff2)-1000)));
            
            halfTrapezoidalMax = trapezoidalMax/2;
            
            [~,newsIndex2] = findpeaks(abs(vOutDiff2(1:numel(vOutDiff2)-1000)), ...
                'MinPeakHeight',halfTrapezoidalMax,'MinPeakDistance', ...
                25E-9/vIn.desc.Ts);
            
            newsIndexBegin = newsIndex2(2);% - round(100*10.^(-9)/vIn.desc.Ts);
            newsIndexEnd = newsIndex2(3);% + round(100*10.^(-9)/vIn.desc.Ts);
            
            %plot(vIn.x,vOut,'-p','MarkerIndices',[newsIndexBegin newsIndexEnd],'MarkerFaceColor','red','MarkerSize',15);
            %orig = vOut;
            %plot(abs(vOutDiff2))
            %baseline = mean(vOut(newsIndex2(1)-round(2E-6/vIn.desc.Ts): ...
            %    newsIndex2(1)));
            
            %vOut = (vOut - baseline);
            %plot(vIn.x,orig,vIn.x,vOut)
            %----------------------------------------------------------------------
            % After pulse baseline correction.
            %----------------------------------------------------------------------
            %{
            if numel(newsIndex) > 4
                secondPulseIndex = newsIndex2(5);
                'shaped pileup'
                
            else
                secondPulseIndex = numel(vOut) - 100;
                
            end
            
            afterPulseBaseline = mean(vOut(round(newsIndex2(4)+ ...
                (100E-9/vIn.desc.Ts)):secondPulseIndex));
            
            if afterPulseBaseline >= 0.01
                trapezoidalPlateauAverage = 0;
                vOut = vInP.y;
                'after'
                return;
                
            else
            %}
                trapezoidalPlateauAverage = mean(vOut(newsIndexBegin: ...
                    newsIndexEnd));
                
            %end
            
        catch
            trapezoidalPlateauAverage = 0;
            
        end
        
        %{
        if trapezoidalPlateauAverage >= 5.93 && trapezoidalPlateauAverage <= 5.964
            figure
            plot(vIn.y);
            title('Higher energy side');
            
        elseif trapezoidalPlateauAverage >= 5.90 && trapezoidalPlateauAverage < 5.93
            figure
            plot(vIn.y);
            title('Lower energy side');
            
        end
        %}
        
        vOut = vOut';
        
        %{
        pFit = polyfit(vIn.x(newsIndexBegin: ...
            newsIndexEnd),vOut(newsIndexBegin:newsIndexEnd),1);
        
        p = polyval(pFit,vIn.x(newsIndexEnd+1000));
        
        %pFit
        
        
        figure
        
        plot(vIn.x(newsIndexBegin: ...
                newsIndexEnd),p,vIn.x(newsIndexBegin: ...
                newsIndexEnd),vOut(newsIndexBegin:newsIndexEnd))
        %}
        
        %movingSum = movsum(difference,L);
        %assignin('base','movingSum',movingSum');
        %assignin('base','difference',difference);
        %vOut = movingSum';
        %diffGauss = vertcat(zeros(2*L,1),diffGauss);
        
        %}
        %----------------------------------------------------------------------
        % Multi-Site Event Cancellation
        %----------------------------------------------------------------------
        %{
    baseline = mean(diffGauss(100:1000));
    [~,deltaIndex] = max(diffGauss(100:numel(diffGauss)-100));
    deltaIndex = deltaIndex + 100;
    sumLeft = 0;
    sumRight = 0;
    foundBaselineLeft = false;
    foundBaselineRight = false;
    i = 1;
    
    %figure;
    %plot(diffGauss);

    while i < (300/vIn.desc.Ts) && deltaIndex - i > 0 && deltaIndex + i < numel(diffGauss)
        
        if diffGauss(deltaIndex-i) > baseline && foundBaselineLeft == false
            sumLeft = sumLeft + diffGauss(deltaIndex-i);
            diffGauss(deltaIndex-i) = 0;
            
        elseif diffGauss(deltaIndex-i) <= baseline
            foundBaselineLeft = true;
        end
        
        if diffGauss(deltaIndex+i) > baseline && foundBaselineRight == false
            sumRight = sumRight + diffGauss(deltaIndex+i);
            diffGauss(deltaIndex+i) = 0;
            
        elseif diffGauss(deltaIndex+i) <= baseline
            foundBaselineRight = true;
        end
        
        i = i + 1;
    end
    
    sumPeak = sumLeft + sumRight;
    diffGauss(deltaIndex) = diffGauss(deltaIndex) + sumPeak;
    
    %figure;
    %plot(diffGauss);
    
    if deltaIndex + L > numel(diffGauss)
        maxSpectrum = numel(diffGauss - 10);
        
    else
        maxSpectrum = deltaIndex + L;
        
    end
    
    if deltaIndex - L <= 0
        minSpectrum = 10;
        
    else
        minSpectrum = deltaIndex - L;
        
    end
        %}
        %{
        minSpectrum = 10;
        maxSpectrum = numel(diffGauss - 10);
        %}
        %{
        for i=1:1
            diffGauss = smooth(diffGauss(minSpectrum:maxSpectrum),L,'moving');
            
        end
        %}
        %figure;
        %plot(diffGauss);
        %{
        for i=1:1
            diffGauss = smooth(diffGauss(10:numel(diffGauss)-10),L-G,'moving');
            
        end
        %}
        %figure;
        %plot(diffGauss);
        %{
    for i=1:2
        diffGauss = smooth(diffGauss(10:numel(diffGauss)-10),L/2,'moving');
        
    end
        %}
        %diffGauss = diffGauss - mean(diffGauss(500:1000));
        
        %vOut = difference*10.^(4);
        
        %trapezoidalPlateauAverage = max(vOut);
        %assignin('base','trapezoidalPlateauAverage',trapezoidalPlateauAverage);
        %save(['J:\trapezoidalSaves\trapezoidal1.txt'], 'vOut','-ascii','-double','-tabs');
        %figure;
        %plot(vOut);
        
        %}
        
        
    elseif shapingTypeIn == 5
        %{
        for i=1:numel(vInP.y)
            if i < L
                cuspFilter(i) = sinh(vIn.x(i));
            elseif i >= L && i < L + G
                cuspFilter(i) = sinh(vIn.x(L));
            elseif i >= L + G && i < 2*L + G
                cuspFilter(i) = sinh(vIn.x(2*L+G-i));
            end
            
        end
        
        %figure
        %plot(cuspFilter)
        assignin('base','cuspFilter',cuspFilter);
        
        vOut = conv(vInP.y,cuspFilter)*10000;
        
        trapezoidalPlateauAverage = max(vOut);
        %}
        %plot(horzcat(pulse',pulse2'));
        %{
        s = evalin('base','s;');
        
        sum1 = 0;
        sum2 = 0;
        a = 0;
        b = 0;
        
        for it=1:30
            for n=1:numel(pulse-1)
                
                sum1 = sum1 + pulse(n)*s(n) - b*s(n);
                sum2 = sum2 + s(n)^2;
                
            end
            
            a = sum1/sum2;
            
            sum1 = 0;
            sum2 = 0;
            
            for n=1:numel(pulse-1)
                
                sum1 = sum1 + pulse(n) - a*s(n);
                
            end
            
            b = sum1/numel(pulse);
            
            sum1 = 0;
            sum2 = 0;
            
        end
        
        sum1 = 0;
        
        x = a*s + b;
        
        %baseline = mean(x(250:300))
        
        err = immse(x,pulse);
        %}
        
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
            vIn.y = vIn.y;% - baseline;
            %[trapezoidalPlateauAverage1,vOut,noGood] = MWD(vIn,TFRFIn,TFTopIn);
            [vIn.y,err,a] = MLSmoothing(vIn.y,Ts);
            err
            
            if err < 2.5E-3
                
                %[trapezoidalPlateauAverage2,vOut,noGood] = MWD(vIn,TFRFIn,TFTopIn);
                
                %trapezoidalPlateauAverage = (trapezoidalPlateauAverage1 + ...
                %    trapezoidalPlateauAverage2 + a*150)/3;
                
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
        %}
        %{
        if a*100 > 2.5 && err < 2.5E-6
            
            vOut = x;
            
            trapezoidalPlateauAverage = a*100 - mean(vIn.y(1000:10000))*100;
            
        elseif a*100 <= 2.5 && err < 1.5E-6
            
            vOut = x;
            
            trapezoidalPlateauAverage = a*100 - mean(vIn.y(1000:10000))*100;
            
        else
            
            vOut = [];
            
            trapezoidalPlateauAverage = 0;
            
        end
        %}
        % -------------------------------------------------------------------
        % Fitting function
        % -------------------------------------------------------------------
        
    elseif shapingTypeIn == 7
        
        xAxis = linspace(1,numel(vIn.y),numel(vIn.y))';
        %vIn.y = smooth(vIn.y,5001,'lowess');
        %plot(vIn.y)
        v = vIn.y(round(10E-6/Ts):round(24E-6/Ts));
        
        [fitResult, gof] = fitting3(xAxis(round(10E-6/Ts):round(24E-6/Ts)),v);
        gof.adjrsquare
        %baseline = mean(vIn.y(10:5000));
        
        double a;
        a = double(vpa(coeffvalues(fitResult)));
        assignin('base','a',a);
        y = feval(fitResult,xAxis);
        amp = (2*a(1))*100;
        %vIns = evalin('base','vIn;');
        
        if abs(gof.adjrsquare) >= 0.998
            
            trapezoidalPlateauAverage = amp;
            vOut = [v*10.^3,feval(fitResult,xAxis(round(10E-6/Ts):round(24E-6/Ts)))*10.^3];
            
            % - baseline)*10.^3;
            %spectrum = fftshift(fft(vIn.y))/N;
            %{
            assignin('base','vIns',v);
            %assignin('base','vInsF',spectrum);
            assignin('base','a',a');
            evalin('base','vIn = horzcat(vIn,vIns);');
            evalin('base','fits = horzcat(fits,a);');
            %}
            
        %elseif abs(gof.adjrsquare) >= 0.995 && amp <= 2.5
            
        %    trapezoidalPlateauAverage = amp;
        %    vOut = [v*10.^3,feval(fitResult,xAxis(round(10E-6/Ts):round(24E-6/Ts)))*10.^3];
            %{
            assignin('base','vIns',v);
            %assignin('base','vInsF',spectrum);
            assignin('base','a',a');
            evalin('base','vIn = horzcat(vIn,vIns);');
            evalin('base','fits = horzcat(fits,a);');
            %}
            %{
        if trapezoidalPlateauAverage <= 13.8 && trapezoidalPlateauAverage >= 13.6
            
            'a'
            vpa(a(1))
            'TPA'
            trapezoidalPlateauAverage
            % Plot fit with data.
            figure( 'Name', 'fitsave copy 2' );
            h = plot( fitResult, xAxis, vIn.y );
            legend( h, 'unshaped vs. t', 'fitsave copy 2', 'Location', 'NorthEast' );
            % Label axes
            xlabel t
            ylabel unshaped
            grid on
            
        end
            %}
            %{
            FFTIn(3,1) = 10*10.^(6);
            FFTIn(3,2) = 10;
            filtered = bpfilter(vIn,3,FFTIn);
            assignin('base','filtered',filtered);
            %}
            %assignin('base','TPA',trapezoidalPlateauAverage);
            %evalin('base','targetValues = horzcat(targetValues,TPA);');
            
        else
            
            vOut = [];
            trapezoidalPlateauAverage = 0;
            
        end
        
        % Parabola fitting for shaped pulse.
        
    elseif shapingTypeIn == 9
        
        [~,pulseMaxIndex] = max(vIn.y);
        
        pulseRangeMin = pulseMaxIndex - L;%round((0.5*10.^(-6))/(vIn.desc.Ts));
        pulseRangeMax = pulseMaxIndex + L;%round((0.5*10.^(-6))/(vIn.desc.Ts));
        
        if pulseRangeMax < numel(vIn.y) && pulseRangeMin > 0
            
            lastwarn('')
            
            [p,S,mu] = polyfit(vIn.x(pulseRangeMin:pulseRangeMax), ...
                vIn.y(pulseRangeMin:pulseRangeMax),2);
            assignin('base','p',p);
            [warnmsg,~] = lastwarn;
            
            [vOut,delta] = polyval(p,vIn.x(pulseRangeMin:pulseRangeMax),S,mu);
            %xValues = vIn.x(pulseRangeMin:pulseRangeMax);
            vOut = [vIn.y(pulseRangeMin:pulseRangeMax),vOut]*100; %p(1)*xValues.^2 + p(2)*xValues + p(3)];
            
            RSS = sum(delta.^2);
            
            baseline = mean(vIn.y(100:100+round(10E-6/Ts)));
            
            if isempty(warnmsg) && RSS <= 2.9*10.^(-4)
                values = p;%coeffvalues(parabolaFit);
                
                trapezoidalPlateauAverage = ((values(3) - values(2)^2/(4*values(1))))*100;% - baseline)*100;
                
            else
                trapezoidalPlateauAverage = 0;
                
            end
            
        else
            trapezoidalPlateauAverage = 0;
            vOut = [];
            
        end
        
        %assignin('base','TPA',trapezoidalPlateauAverage);
        %evalin('base','targetValues = horzcat(targetValues,TPA);');
        
        %{
    figure( 'Name', 'fitsave copy 2' );
    h = plot( p, vIn.x, vIn.y );
    legend( h, 'unshaped vs. t', 'fitsave copy 2', 'Location', 'NorthEast' );
    % Label axes
    xlabel t
    ylabel unshaped
    grid on
        %}
        
    elseif shapingTypeIn == 10
        %{
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
        %}
        %baseline = mean(vIn.y(round((offset-999E-9)*Fs):round((offset-500E-9)*Fs)));
        
        %v = vIn.y(round(10E-6/Ts)+240:round(30E-6/Ts)-1084) - baseline;
        %size(v)
        %v = v - mean(v(100:250));
        %assignin('base','vIns',vIn.y(round(10E-6/Ts):round(27E-6/Ts)-771));
        %assignin('base','vIns2',v);
        %assignin('base','vIns3',vIn.y(round(10E-6/Ts)+399:round(27E-6/Ts)-771));
        %aL = evalin('base','net(vIns2);');
        %assignin('base','aL',aL);
        %aP = evalin('base','net(vIns2);');
        %aS = evalin('base','net(vIns2);');
        %assignin('base','a',a);
        %aL = evalin('base','net(a);');
        %a(2,:) = evalin('base','net5(vIns);');
        %a(3,:) = evalin('base','net5(vIns);');
        %a(4,:) = evalin('base','net4(vIns);');
        %a2 = a(2,:);
        %a3 = a(3,:);
        %aF = (aL + aP)./2;% + a2 + a3)./3;
        %plot([a2 a2]) 
        %if a >= 3.5
        %plot([v aL])
        %err = immse(aL,v)
        
        %aF = aL; %smooth(aL,100,'moving');
        %aF = smooth(aF,175,'moving');
        %aF = smooth(aF,175,'moving');
        %aF2 = smoothdata(aP,'gaussian');
        %aF3 = smoothdata(aS,'rloess');
        
        %for l=1:5
        %    aF = smooth(aF,50,'moving');
            
        %end
        %aF = smooth(aF,150,'moving');
        
        %if err <= 1E-6
        
        %p = vIn.y(round(16.25E-6/Ts+3):round(35E-6/Ts)) + 70;
        
        %plot(p)
        p = vIn.y(round(10E-6/Ts):round(35E-6/Ts));
        size(p)
        %p = p(400:2000);
        %net = evalin('base','net4');
        
        pM = evalin('base','ps1Mean');
        pS = evalin('base','ps1Std');
        %out = (predict(net, , 'ExecutionEnvironment', 'gpu'));
        
        %plot([(p - pM)/pS])
        %baseline = 0;%mean(out(1:50));
        %out = out - baseline;
        %}
        p = p(251:2000);
        net = evalin('base','net4');
        %spec = spectrogram(p, 100, 'yaxis');
        %rispec = horzcat(real(spec(1:70,:)), imag(spec(1:70,:)));
        
        out = predict(net, (p - pM)/pS, 'ExecutionEnvironment', 'gpu');
        
        out = out - mean(out(1:200));
        
        plot([(p - pM)/pS, out])
        trapezoidalPlateauAverage = (max(smooth(out,25,'moving')))  
        %trapezoidalPlateauAverage = (max(aF(1:numel(aF)-10))*100 + max(aF2(1:numel(aF2)-10))*100 + ...
        %    max(aF3(1:numel(aF3)-10))*100)/3;% - mean(vIn.y(round(10E-6/Ts):round(10E-6/Ts)+400))*100;% - mean(a(100:400))*100;
        vOut = out;
        noGood = 0;
        
        if trapezoidalPlateauAverage <= 0 && false
            trapezoidalPlateauAverage = 0;
            vOut = [];
            noGood = 1;
            %err
            
        end
        %plot([aF vIn.y(round(10E-6/Ts)+299:round(24E-6/Ts)-559)]);
        %if false && trapezoidalPlateauAverage <= 1.4
        %plot([aL aF v]);
        %vOut = [];
        %trapezoidalPlateauAverage = 0;
        %end
        %else
        %   trapezoidalPlateauAverage = 0;
        %   vOut = [];
        %end
        
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
        %{
        FFTIn(3,1) = 10*10.^(6);
        FFTIn(3,2) = 10;
        filtered = bpfilter(vIn,3,FFTIn);
        assignin('base','filtered',filtered);
        trapezoidalPlateauAverage = evalin('base','(NNet{20}(filtered)+NNet{11}(filtered))/2;');
        vOut = [];
        %}
    end
    
else
    trapezoidalPlateauAverage = 0;
    vOut = [];
    
end
%{
catch
    trapezoidalPlateauAverage = 0;
    vOut = [];
    noGood = 1;
    'trapezoidal filter error'
    
end
%}
end