function [binsRiseTime,riseTimeHistogram] = RiseTimeHistogram(numSamples,... 
    selectPath,FFTIn,numChannel,BPF)
tic
c = numChannel;
FFTIn(numChannel,1) = FFTIn(numChannel,1)*10.^(6);
FFTIn(numChannel,2) = FFTIn(numChannel,2)*10.^(6);
d = [];

for i=1:numSamples
    stop1 = 0;
    stop2 = 0;
    count1 = 0;
    count2 = 0;
    channelNumber = cell(2,1);
    pulse = cell(1,2);
    foundRise1 = 1;
    foundRise2 = 1;
    pulseIny = 0;
    pulseInx = 0;
    
    if (true)
        channelNumber{c} = ['C',stringconversion(c),'1'];
        %appendage = num2str(i-1,'%05.f');
        
        if i < 10
            appendage = ['0000',stringconversion(i)];
        elseif i >= 10 && i < 100
            appendage = ['000',stringconversion(i)];
        elseif i >= 100 && i < 1000
            appendage = ['00',stringconversion(i)];
        elseif i >= 1000 && i < 10000
            appendage = ['0',stringconversion(i)];
        else
            appendage = stringconversion(i);
        end
    end
    
    if (false)
        channelNumber{c} = ['C',stringconversion(c),'XX'];
        %appendage = num2str(i-1,'%05.f');
        
        if i < 10
            appendage = ['000000000',stringconversion(i)];
        elseif i >= 10 && i < 100
            appendage = ['00000000',stringconversion(i)];
        elseif i >= 100 && i < 1000
            appendage = ['0000000',stringconversion(i)];
        elseif i >= 1000 && i < 10000
            appendage = ['000000',stringconversion(i)];
        elseif i >= 10000 && i < 100000
            appendage = ['00000',stringconversion(i)];
        elseif i >= 100000 && i < 1000000
            appendage = ['0000',stringconversion(i)];
        elseif i >= 1000000 && i < 10000000
            appendage = ['000',stringconversion(i)];
        elseif i >= 10000000 && i < 100000000
            appendage = ['00',stringconversion(i)];
        elseif i >= 10000000 && i < 100000000
            appendage = ['0',stringconversion(i)];
        else
            appendage = stringconversion(i);
        end
        
    end
    
    fileName = [selectPath,'\',channelNumber{c},appendage,'.trc'];
    
    pulse{c}(i) = waveform.ReadLeCroyBinaryWaveform(fileName);
    
    if BPF == 1
        [pulse{c}(i).y,d] = bpfilter(pulse{c}(i),c,FFTIn,d);
        
    end
    
    pulseIn.y = pulse{c}(i).y;
    pulseIn.x = pulse{c}(i).x;

    [VMin,VMinIndex] = max(abs(pulseIny(50:numel(pulseIny)-50)));
    VMinIndex = VMinIndex + 50;
    riseTimeLower = (0.1*VMin);
    riseTimeHigher = (0.9*VMin);
    
    % Manual rise time calculation.
    %{
    while (stop1 == 0 || stop2 == 0) && VMinIndex > count2 && VMinIndex > count1 ...
            && count2 < 3000 && count1 < 3000
        
        % Going down peak edge by 'count' until less than VMinFraction.
        if stop1 == 0 && abs(pulseIny(VMinIndex - count1)) <= abs(riseTimeLower)
            foundRise1 = (VMinIndex - count1);
            stop1 = 1;
        elseif stop1 == 0
            count1 = (count1 + 1);
        end
        
        % Going down peak edge by 'count' until less than VMinFraction.
        if  stop2 == 0 && abs(pulseIny(VMinIndex - count2)) <= abs(riseTimeHigher)
            foundRise2 = (VMinIndex - count2);
            stop2 = 1;
        elseif stop2 == 0
            count2 = (count2 + 1);
        end
    end
    %}
    %paes(timingIn,TypeIn,c,multiStopIn,multiStopConditionsIn,numPeaks,pulseIn,RTFLIn,RTFUpIn)
    
    [~,~,~,riseTimeCalculate] = paes(0,3,c,'0',1,1,pulse{c}(i),0,1000);
    
    riseTime(i) = riseTimeCalculate; %pulseInx(foundRise2) - pulseInx(foundRise1);
    
    %{
    riseInputY = pulseIny(500:numel(pulseIny)-500);
    riseInputX = pulseInx(500:numel(pulseIny)-500);
    riseTimeCalculate(i) = stepinfo(riseInputY,riseInputX);
    riseTime(i) = riseTimeCalculate(i).RiseTime;
    %}
end
assignin('base','riseTime',riseTime);
binsRiseTime = linspace(0,20E-9,100);

riseTimeHistogram = hist(riseTime,binsRiseTime);
assignin('base','riseTimeHistogram',riseTimeHistogram);
riseTimeHistogram(100) = 0;
riseTimeHistogram(1) = 0;
binsRiseTime = binsRiseTime';
riseTimeHistogram = riseTimeHistogram';
toc