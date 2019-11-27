function [output,value,maximum] = neuralNetworkFiltering(input,timeIndex)

try
    smoothed = input.y; %smooth(input.y,251,'lowess');
    
    value = smoothed(timeIndex+round((25E-9/input.desc.Ts)));
    
    smoothReversed = smoothed;
    
    %figure
    %plot(smoothReversed)
    
    [maximum,maximumIndex] = max(smoothReversed(100:numel(smoothReversed)-100));
    
    if timeIndex < maximumIndex
        [~,beforeMax] = max(input.y(100:numel(input.y)-100));
        beforeMax = beforeMax + 100;
        maximumIndex = maximumIndex + 100;
        %maximumIndex+round(0.4E-6/input.desc.Ts)
        
        try
            output = smoothReversed(maximumIndex-round(0.5E-6/input.desc.Ts):maximumIndex+round(0.2E-6/input.desc.Ts));
        catch
            smoothReversed = padarray(smoothReversed,round(0.5E-6/input.desc.Ts),'pre');
            output = smoothReversed(maximumIndex-round(0.5E-6/input.desc.Ts):maximumIndex+round(0.2E-6/input.desc.Ts));
        end
        %figure
        %plot(output)
        
        valueIndex = timeIndex+round((25E-9/input.desc.Ts));
        
        %plot(smoothed,'-p','MarkerIndices',valueIndex,'MarkerFaceColor','red','MarkerSize',15);
    else
        output = [];
        value = 0;
        maximum = 0;
        
    end
catch
    output = [];
    value = 0;
    maximum = 0;
    
end

end