numBins = 512;
cutoff = 810E-9;
timingWindow = 1E-6;

eVse = zeros(numBins,numBins);

binsSize = (timingWindow)/numBins;

for s=1:length(totalMultiToFs(1,:))
    
    whichBinToF(1,s) = ceil(totalMultiToFs(1,s)/binsSize);
    whichBinToF(2,s) = ceil(totalMultiToFs(2,s)/binsSize);
    
    %{
    if s > 1
    'bin ratio'
    (whichBinToF(2,s) - whichBinToF(2,s-1))/(whichBinToF(1,s) - whichBinToF(1,s-1))
    end
    %'bin2'
    %whichBinToF(2,s)
    %}
    if whichBinToF(1,s) > 0 && whichBinToF(1,s) <= numBins && ...
            whichBinToF(2,s) > 0 && whichBinToF(2,s) <= numBins && ...
            whichBinToF(1,s) < whichBinToF(2,s) %&& ...
            %totalMultiToFs(2,s) <= cutoff
        
        eVse(whichBinToF(1,s),whichBinToF(2,s)) = eVse(whichBinToF(1,s),whichBinToF(2,s)) + 1;
        
    end
    
end