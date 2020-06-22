function [hToF,hGe,geMaxVsToF,gVsg,eVse] = histo(numberOfSamples,multiStop,numBins,ToFChannel, ...
    gammaWindow,gammaLower,timingWindow,geOutMax,timeOfFlight,energyCoinCha,inputChannel,GGCoinCha, ...
    GGCoinOn,eeCoinOn,multiStopConditionsIn,secondECutoff)

if GGCoinOn == false && eeCoinOn == false
    
    % --------------------------------------------------------------------------
    % Histogram creation function for time-of-flight, gamma, and 2D spectra.
    % --------------------------------------------------------------------------
    
    whichBinGe = zeros(1,numberOfSamples);
    
    % Creating histogram for the time of flight spectrum.
    if isempty(timeOfFlight) == 0 && ToFChannel ~= 5
        if multiStop == '1'
            whichBinTimeOfFlight = cell(1,numberOfSamples);
            
        else
            whichBinTimeOfFlight = zeros(1,numberOfSamples);
            
        end
        
        if timingWindow > 1
            ToFInitialTime = 0;
            
        else
            ToFInitialTime = -250E-9;
            
        end
        
        %bins = linspace(ToFInitialTime,timingWindow,numBins(ToFChannel))';
        binSize = (timingWindow - ToFInitialTime)/numBins(ToFChannel);
        
        if multiStop == '1'
            timeIn1 = 1;
            assignin('base','timeOfFlight',timeOfFlight);
            for q=1:numberOfSamples
                %size(timeOfFlight{q})
                %size(timeIn1)
                timeIn1 = horzcat(timeIn1,timeOfFlight{q});
                
            end
            
            hToF = histcounts(timeIn1,'NumBins',numBins(ToFChannel),'BinLimits', ...
                [ToFInitialTime timingWindow])';
            %hToF = hToFHist.Values';
            
        else
            %'IN HISTOGRAM'
            %timeOfFlight
            %pause(2)
            hToF = histcounts(timeOfFlight,'NumBins',numBins(ToFChannel),'BinLimits', ...
                [ToFInitialTime timingWindow])';
            %assignin('base','hToF',hToF);
            %hToF = hToFHist.Values';
            
        end
        
    else
        hToF = [];
        
    end
    
    % Creating histogram for the gamma spectrum.
    if isempty(geOutMax) == 0 && inputChannel ~= ToFChannel
        VoutMax2 = geOutMax;
        
        bins3 = linspace(gammaLower,gammaWindow,numBins(inputChannel))';
        bin3Size = ((gammaWindow - gammaLower)/numBins(inputChannel));
        
        hGe = histcounts(VoutMax2,'NumBins',numBins(inputChannel),'BinLimits',[gammaLower gammaWindow])';
        %hGe = hGeHist.Values';
        
    else
        hGe = [];
        
    end
    
    % 3-dimensional plot of number of counts vs gamma energy and time of flight
    if numel(timeOfFlight) > 0 && numel(geOutMax) > 0 && multiStop == '1' && ...
            energyCoinCha == inputChannel && ToFChannel ~= 5
        
        geMaxVsToF = zeros(numBins(energyCoinCha),numBins(ToFChannel));
        
        for s=1:numberOfSamples
            whichBinGe(s) = floor(VoutMax2(s)/bin3Size);
            
            if multiStopConditionsIn ~= 6
                
                for q=1:numel(timeOfFlight{s})
                    whichBinTimeOfFlight{s}(q) = floor((timeOfFlight{s}(q)+250*10.^(-9))/binSize);
                    
                    if whichBinGe(s) > 0 && whichBinTimeOfFlight{s}(q) > 0 && ...
                            whichBinTimeOfFlight{s}(q) <= numBins(ToFChannel) && ...
                            whichBinGe(s) <= numBins(energyCoinCha) && ...
                            multiStopConditionsIn ~= 6
                        
                        geMaxVsToF(whichBinGe(s),whichBinTimeOfFlight{s}(q)) =  ...
                            geMaxVsToF(whichBinGe(s),whichBinTimeOfFlight{s}(q)) + 1;
                        
                    end
                    
                end
            end
        end
        
    elseif numel(timeOfFlight) > 0 && numel(geOutMax) > 0 && multiStop == '0' && ...
            energyCoinCha == inputChannel
        
        geMaxVsToF = zeros(numBins(energyCoinCha),numBins(ToFChannel));
        
        for s=1:numberOfSamples
            
            whichBinGe(s) = floor(VoutMax2(s)/bin3Size);
            whichBinTimeOfFlight(s) = floor((timeOfFlight(s)+250*10.^(-9))/binSize);
            
            if whichBinGe(s) > 0 && whichBinTimeOfFlight(s) > 0 && whichBinTimeOfFlight(s) <= numBins(ToFChannel) ...
                    && whichBinGe(s) <= numBins(energyCoinCha)
                
                geMaxVsToF(whichBinGe(s),whichBinTimeOfFlight(s)) = geMaxVsToF(whichBinGe(s),whichBinTimeOfFlight(s)) + 1;
                
            end
        end
        
    else
        geMaxVsToF = [];

%{
figure
surf(electronVsElectron);
shading interp;
colormap(jet);
rotate3d on;
xlim([1 4000]);
ylim([1 4000]);
zlim([1 50]);
%set(plo,'LineStyle','none');
%mesh(geMaxVsTimeOfFlight);
view([45 0]);
%daspect([1 1 1]);
%camzoom(1.4);

%imagesc(flipud(geMaxVsTimeOfFlight));
%axis tight;
%set(gca,'YDir','normal');

toc
        %}
    end
    
    gVsg = [];
    eVse = [];
    
elseif GGCoinOn == true && eeCoinOn == false
    
    gVsg = zeros(numBins(GGCoinCha(1)),numBins(GGCoinCha(2)));
    
    for g=1:2
        bins(g,:) = linspace(gammaLower(GGCoinCha(g)),gammaWindow(GGCoinCha(g)),numBins(GGCoinCha(g)))';
        binsSize(g) = ((gammaWindow(GGCoinCha(g)) - gammaLower(GGCoinCha(g)))/numBins(GGCoinCha(g)));
        
    end
    
    for s=1:numberOfSamples
        assignin('base','geOutMax',geOutMax);
        whichBinGe(1,s) = floor(geOutMax{GGCoinCha(1)}(s)/binsSize(1));
        whichBinGe(2,s) = floor(geOutMax{GGCoinCha(2)}(s)/binsSize(2));
        
        if whichBinGe(1,s) > 0 && whichBinGe(GGCoinCha(1),s) <= numBins(GGCoinCha(1)) && ...
                whichBinGe(2,s) > 0 && whichBinGe(GGCoinCha(2),s) <= numBins(GGCoinCha(2))
            
            gVsg(whichBinGe(1,s),whichBinGe(2,s)) = gVsg(whichBinGe(1,s),whichBinGe(2,s)) + 1;
            
        end
        
    end
    
    hToF = [];
    geMaxVsToF = [];
    hGe = [];
    eVse = [];
    
elseif GGCoinOn == false && eeCoinOn == true
    VoutMax2 = geOutMax;
    gVsg = [];
    'CUTOFF'
    cutoff = double(secondECutoff)*10.^(-9);
    
    if energyCoinCha == inputChannel
        geMaxVsToF = zeros(numBins(energyCoinCha),numBins(ToFChannel),numBins(ToFChannel),'uint16');
        
    end
    
    if numel(timeOfFlight) > 0 && numel(geOutMax) > 0 && multiStop == '1' && ...
            energyCoinCha == inputChannel && ToFChannel ~= 5
        
        if timingWindow > 1
            ToFInitialTime = 0;
            
        else
            ToFInitialTime = -250E-9;
            
        end
        
        bin3Size = ((gammaWindow - gammaLower)/numBins(inputChannel));
        binSize = (timingWindow - ToFInitialTime)/numBins(ToFChannel);
        
        for s=1:numberOfSamples
            whichBinGe(s) = ceil((VoutMax2(s) - gammaLower)/bin3Size);
            
            if multiStopConditionsIn == 6 && numel(timeOfFlight{s}) == 2 && ...
                    ~isnan(timeOfFlight{s}(1)) && ~isnan(timeOfFlight{s}(2))
                'ASDFLASGKAFGKDFGASGASFG'
                cutoff
                whichBinTimeOfFlight{s}(1) = ceil((timeOfFlight{s}(1)+250*10.^(-9))/binSize);
                whichBinTimeOfFlight{s}(2) = ceil((timeOfFlight{s}(2)+250*10.^(-9))/binSize);
                
                if whichBinGe(s) > 0 && whichBinTimeOfFlight{s}(1) > 0 && ...
                        whichBinTimeOfFlight{s}(1) <= numBins(ToFChannel) && ...
                        whichBinGe(s) > 0 && whichBinTimeOfFlight{s}(2) > 0 && ...
                        whichBinTimeOfFlight{s}(2) <= numBins(ToFChannel) && ...
                        whichBinGe(s) <= numBins(energyCoinCha) && ...
                        timeOfFlight{s}(2) <= cutoff
                        
                    
                    geMaxVsToF(whichBinGe(s),whichBinTimeOfFlight{s}(1),whichBinTimeOfFlight{s}(2)) =  ...
                        geMaxVsToF(whichBinGe(s),whichBinTimeOfFlight{s}(1),whichBinTimeOfFlight{s}(2)) + 1;
                    
                end
                
        %assignin('base','geMaxVsToF',geMaxVsToF);
            end
        end
        
        eVse = [];
        
    end
    
    if energyCoinCha ~= inputChannel
        
        eVse = zeros(numBins(ToFChannel),numBins(ToFChannel));
        
        binsSize = (timingWindow)/numBins(ToFChannel);
        
        for s=1:numberOfSamples
            
            whichBinToF(1,s) = ceil(timeOfFlight{1}(s)/binsSize);
            whichBinToF(2,s) = ceil(timeOfFlight{2}(s)/binsSize);
            
            if whichBinToF(1,s) > 0 && whichBinToF(1,s) <= numBins(ToFChannel) && ...
                    whichBinToF(2,s) > 0 && whichBinToF(2,s) <= numBins(ToFChannel)
                
                eVse(whichBinToF(1,s),whichBinToF(2,s)) = eVse(whichBinToF(1,s),whichBinToF(2,s)) + 1;
                
            end
            
        end
        
        geMaxVsToF = [];
        
    end
    
    hToF = [];
    hGe = [];
    gVsg = [];
    
end

end
