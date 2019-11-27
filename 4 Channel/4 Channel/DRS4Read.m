function [amp,shapedPulse] = DRS4Read(selectPath,c)

fid = selectPath;

numSamples = 1;

findingTCh1 = fseek(fid,16,'bof');

readingTCh1 = fread(fid,1024,'single');

%findingTCh2 = fseek(fid,4116,'bof');

%readingTCh2 = fread(fid,1024,'single');

%readingVCh1 = cell([1,numSamples]);

%readingVCh2 = cell([1,numSamples]);

%readingtCell = zeros(1,numSamples);

%rangeCenter = zeros(1,numSamples);

readingInitialTime = 0;

readingFinalTime = 0;

findingInitialTime = fseek(fid,8224,'bof');

readingInitialTime = fread(fid,4,'short');

findingFinalTime = fseek(fid,(numSamples-1)*4144 + 8224,'bof');

readingFinalTime = fread(fid,4,'short');

disp('Initial time: ');

disp(readingInitialTime);

disp('Final time: ');

disp(readingFinalTime);
c
findingVCh1 = fseek(fid,(c-1)*4176+4152,'bof');

readingVCh1 = fread(fid,1024,'short');

%findingVCh2 = fseek(fid,(c-1)*4144+10308,'bof');

%readingVCh2 = fread(fid,1024,'short');

findingtCell = fseek(fid,(c-1)*6184+4146,'bof');

readingtCell = fread(fid,1,'short');

findingRC = fseek(fid,(c-1)*4176+4130,'bof');

rangeCenter = fread(fid,1,'short')
%{
for i=1:1024
    if readingVCh1(i) < 0
        readingVCh1(i) = ((rangeCenter + 0.95)*65535);
    end
    %{
        if readingVCh2{c}(i) < 0
            readingVCh2{c}(i) = ((rangeCenter(c) + 0.5)*65535);
        end
    %}
end
%}
%}
% Converting to time.
%{
timeConversionCh1 = cell([1,numSamples]);

for c=1:numSamples
    timeConversionCh1Temp = zeros(1,1024);
    'time conversion'
    for i=1:1024
        for j=1:i-1
            timeConversionCh1Temp(i) = (timeConversionCh1Temp(i) + readingTCh1(mod(j+readingtCell(c),1024)+1));
        end
    end
    timeConversionCh1 = timeConversionCh1Temp;
end
%}
%{
timeConversionCh2 = cell([1,numSamples]);

for c=1:numSamples
    timeConversionCh2Temp = zeros(1,1024);
    for i=1:1024
        for j=1:i-1
            timeConversionCh2Temp(i) = (timeConversionCh2Temp(i) + readingTCh2(mod(j+readingtCell(c),1024)+1));
        end
    end
    timeConversionCh2{c} = timeConversionCh2Temp;
end
%}
% Converting to voltage.

'voltage conversion'
for c=1:numSamples
    for i=1:1024
        voltageConversionCh1(i) = ((readingVCh1(i)/65535) - 0.05);
    end
end
%{
voltageConversionCh2 = cell([1,numSamples]);

for c=1:numSamples
    for i=1:1024
        voltageConversionCh2{c}(i) = ((readingVCh2{c}(i)/65535) - (rangeCenter(c) + 0.5));
    end
end
%}

% - - - - - - - - - Voltage baseline correction - - - - - - - - -
%{
sum1 = zeros(1,numSamples);
sum2 = zeros(1,numSamples);
count1 = zeros(1,numSamples);
count2 = zeros(1,numSamples);
baselineAverageCh1 = zeros(1,numSamples);
baselineAverageCh2 = zeros(1,numSamples);

for c=1:1
    for i=200:250
        if voltageConversionCh1(i) <= 0
            sum1(c) = (sum1(c) + voltageConversionCh1(i));
            count1(c) = (count1(c) + 1);
        end
    end
    
    baselineAverageCh1(c) = (sum1(c)/count1(c));

end

for c=1:1
    for i=1:1024
        voltageConversionCh1BLCorr(i) = (voltageConversionCh1(i) - baselineAverageCh1(c));
    end
end
%}
waveform = voltageConversionCh1;
%waveform.x = timeConversionCh1;

%{
for c=1:1
    for i=1:200
        if voltageConversionCh2{c}(i) <= 0
            sum2(c) = (sum2(c) + voltageConversionCh2{c}(i));
            count2(c) = (count2(c) + 1);
        end
    end
    
    baselineAverageCh2(c) = (sum2(c)/count2(c));    
    
end

voltageConversionCh2BLCorr = cell([1,numSamples]);

for c=1:numSamples
    for i=1:1024
        voltageConversionCh2BLCorr{c}(i) = (voltageConversionCh2{c}(i) - baselineAverageCh2(c));
    end
end
%}
%{
% Time shifting Channel 2 to align both channels.

timeZeroCh1 = zeros(1,numSamples);
timeZeroCh2 = zeros(1,numSamples);
timeDiff = zeros(1,numSamples);

for c=1:numSamples
    if readingtCell(c) > 1
        timeZeroCh1(c) = timeConversionCh1{c}(mod(1024-readingtCell(c),1024));
        timeZeroCh2(c) = timeConversionCh2{c}(mod(1024-readingtCell(c),1024));
        timeDiff(c) = (timeZeroCh1(c) - timeZeroCh2(c));
    else
        timeZeroCh1(c) = 10;
        timeZeroCh2(c) = 1;
        timeDiff(c) = (timeZeroCh1(c) - timeZeroCh2(c));
    end
end

timeShiftCh2 = cell([1,numSamples]);

for c=1:numSamples
    for i=1:1024
        timeShiftCh2{c}(i) = (timeConversionCh2{c}(i) + timeDiff(c));
    end
end

trueTimeCh1 = cell([1,numSamples]);
trueTimeCh2 = cell([1,numSamples]);

for c=1:numSamples
    trueTimeCh1{c} = timeConversionCh1{c};
    trueTimeCh2{c} = timeShiftCh2{c};
end

% - - - - - - - - - CFD Work - - - - - - - - -

findFrac1 = zeros(1,numSamples);
foundFrac1 = zeros(1,numSamples);
findFrac2 = zeros(1,numSamples);
foundFrac2 = zeros(1,numSamples);
VMinCh1 = zeros(1,numSamples);
VMinCh2 = zeros(1,numSamples);
VMinCh1Index = zeros(1,numSamples);
VMinCh2Index = zeros(1,numSamples);
VMinCh1Time = zeros(1,numSamples);
VMinCh2Time = zeros(1,numSamples);
VMinCh1Fraction = zeros(1,numSamples);
VMinCh2Fraction = zeros(1,numSamples);

for c=1:numSamples
    [VMinCh1(c),VMinCh1Index(c)] = max(abs(voltageConversionCh1BLCorr{c}(100:1024)));
    [VMinCh2(c),VMinCh2Index(c)] = max(abs(voltageConversionCh2BLCorr{c}(100:1024)));

    VMinCh1Time(c) = trueTimeCh1{c}(VMinCh1Index(c));
    VMinCh2Time(c) = trueTimeCh2{c}(VMinCh2Index(c));
    
    % Fraction for CFD.
    
    VMinCh1Fraction(c) = (0.30*VMinCh1(c));
    VMinCh2Fraction(c) = (0.30*VMinCh2(c));

    findFrac1(c) = VMinCh1Index(c);
    findFrac2(c) = VMinCh2Index(c);
end

% Charging from the peak towards lower time to find the closest value to
% the fractional voltage.

timeOfFlight = zeros(1,numSamples);
fractionCh1TimeConversion = zeros(1,numSamples);
fractionCh2TimeConversion = zeros(1,numSamples);
count1 = zeros(1,numSamples);
count2 = zeros(1,numSamples);
counter = zeros(1,numSamples);
    
for c=1:numSamples
    
    % To remove over-pulses, voltages should be less than 200 mV:
    if VMinCh1(c) < 0.2 && VMinCh2(c) < 0.2
    stop = 0;
    
    % Ensuring pulse indices are greater than 1.
    if findFrac1(c) > 1
    while stop == 0 && findFrac1(c) > count1(c)
        if abs(voltageConversionCh1BLCorr{c}(findFrac1(c) - count1(c))) <= abs(VMinCh1Fraction(c))
            foundFrac1(c) = (findFrac1(c) - count1(c));
            stop = 1;
        else
            count1(c) = (count1(c) + 1);
        end
    end
    end

    stop = 0;
    
    if findFrac2(c) > 1
    while stop == 0 && findFrac2(c) > count2(c)
        if abs(voltageConversionCh2BLCorr{c}(findFrac2(c) - count2(c))) <= abs(VMinCh2Fraction(c))
            foundFrac2(c) = (findFrac2(c) - count2(c));
            stop = 1;
        else
            count2(c) = (count2(c) + 1);
        end
    end
    end
    

    
    if findFrac1(c) - count1(c) > 0 && findFrac2(c) - count2(c) > 0 &&  trueTimeCh1{c}(findFrac1(c)) > trueTimeCh2{c}(findFrac2(c))
        
    fractionCh1TimeConversion(c) = trueTimeCh1{c}(findFrac1(c) - count1(c));

    interpTimeCh1{c} = trueTimeCh1{c}(findFrac1(c) - count1(c)):0.01:trueTimeCh1{c}(findFrac1(c) - count1(c) + 1);
    
    interpolationCh1{c} = interp1(trueTimeCh1{c},voltageConversionCh1BLCorr{c},interpTimeCh1{c},'spline');

    [fractionCh1(c), fractionCh1Index(c)] = min(abs(abs(interpolationCh1{c}) - abs(VMinCh1Fraction(c))));
    
    
    fractionCh2TimeConversion(c) = trueTimeCh2{c}(findFrac2(c) - count2(c));
    
    interpTimeCh2{c} = trueTimeCh2{c}(findFrac2(c) - count2(c)):0.01:trueTimeCh2{c}(findFrac2(c) - count2(c) + 1);
    
    interpolationCh2{c} = interp1(trueTimeCh2{c},voltageConversionCh2BLCorr{c},interpTimeCh2{c},'spline');

    [fractionCh2(c), fractionCh2Index(c)] = min(abs(abs(interpolationCh2{c}) - abs(VMinCh2Fraction(c))));
    
    % Finding time difference between fractional values of the peaks.
    
    % timeOfFlight(c) = abs(fractionCh1TimeConversion(c) - fractionCh2TimeConversion(c));
     timeOfFlight(c) = (interpTimeCh1{c}(fractionCh1Index(c)) - interpTimeCh2{c}(fractionCh2Index(c)));
     
    else
        timeOfFlight(c) = 0;
    end
    
    else
        timeOfFlight(c) = 0;
    end

end

bins = zeros(1,2048);

for i=2:2048
    bins(i) = (bins(i-1) + 0.4);
end

h = hist(timeOfFlight,bins);
%}

%waveform.y = waveform.y(200:numel(waveform.y));
%waveform.x = waveform.x(200:numel(waveform.x));

shapedPulse = waveform;%smooth(waveform.y,351,'moving');



amp = max(shapedPulse);

xAxis = linspace(1,numel(shapedPulse),numel(shapedPulse))';
        %vIn.y = smooth(vIn.y,5001,'lowess');
        %plot(vIn.y)
        v = shapedPulse';
        
        [fitResult, gof] = fitting4(xAxis,v);
        gof.adjrsquare
        %baseline = mean(vIn.y(10:5000));
        
        double a;
        a = double(vpa(coeffvalues(fitResult)));
        assignin('base','a',a);
        y = feval(fitResult,xAxis);
        amp = (2*a(1))*100;
        %vIns = evalin('base','vIn;');
            plot([y v])
            
            if gof.adjrsquare < 0.99
                amp = NaN;
                
            end


%{
        [~,pulseMaxIndex] = max(shapedPulse);
        
        pulseRangeMin = pulseMaxIndex - 20;%round((0.5*10.^(-6))/(vIn.desc.Ts));
        pulseRangeMax = pulseMaxIndex + 20;%round((0.5*10.^(-6))/(vIn.desc.Ts));
        x = linspace(0,1024,1024);
        
        if pulseRangeMax < numel(shapedPulse) && pulseRangeMin > 0
            
            lastwarn('')
            
            [p,S,mu] = polyfit(x(pulseRangeMin:pulseRangeMax), ...
                shapedPulse(pulseRangeMin:pulseRangeMax),2);
            assignin('base','p',p);
            [warnmsg,~] = lastwarn;
            
            [vOut,delta] = polyval(p,x(pulseRangeMin:pulseRangeMax),S,mu);
            %xValues = vIn.x(pulseRangeMin:pulseRangeMax);
            vOut = [shapedPulse(pulseRangeMin:pulseRangeMax),vOut]*100; %p(1)*xValues.^2 + p(2)*xValues + p(3)];
            
            RSS = sum(delta.^2);
            
            %baseline = mean(shapedPulse(100:100+round(10E-6/Ts)));
            
            if isempty(warnmsg) && RSS <= 1%2.9*10.^(-4)
                values = p;%coeffvalues(parabolaFit);
                
                amp = ((values(3) - values(2)^2/(4*values(1))))*100;% - baseline)*100;
                
            else
                amp = 0;
                
            end
            
        else
            amp = 0;
            vOut = [];
            
        end
        %}
end