xToF = linspace(-1E-9,1.5E-9,256)';


    observe = [];
    fwhm = [];
    index = [];
    best2 = [];
    amplitude = [];
    x0 = [];
    

for n=1:64
    
    %if exact(n) == 1
        
        for i=1:15
            
            if allNeurons.eletParameters(1,i) ~= 0
                
                allNeurons.neurons{n}(allNeurons.neurons{n}(:,i) == 0,i) = NaN;
                
                hToF = histcounts(allNeurons.neurons{n}(:,i),'NumBins',256,'BinLimits',[-1E-9 1.5E-9])';
                hToF = hToF/max(hToF);
                plot(hToF)
                observe(i,n) = numel(allNeurons.neurons{n}(~isnan(allNeurons.neurons{n}(:,i)),i));
                
                if max(hToF) < sum(hToF)
                    
                    try
                        
                        [~,second] = max(hToF);
                        
                        gaussFit = fit(xToF,hToF,'gauss1');
                        
                        cValues = coeffvalues(gaussFit);
                        
                        
                        fwhm(i,n) = (2*sqrt(log(2))*cValues(3))/sqrt(2);
                        amplitude(i,n) = cValues(1);
                        x0(i,n) = cValues(2);
                        
                        if fwhm(i,n) == 0
                            fwhm(i,n) = 1;
                            
                        end
                        
                        %{
                        nVector = allNeurons.neurons{n}(isnan(allNeurons.neurons{n}(:,i)) == 0 | ~any(allNeurons.neurons{n}(:,i)),i);
                        r = fitdist(nVector,'Normal');
                        fwhm(i,n) = r.sigma;
                        x0(i,n) = r.mu;
                        %}
                        %{
                        nVector = allNeurons.neurons{n}(isnan(allNeurons.neurons{n}(:,i)) == 0 | ~any(allNeurons.neurons{n}(:,i)),i);
                        [mu,sigma] = normfit(nVector);
                        fwhm(i,n) = sigma;
                        x0(i,n) = mu;
                        %}
                        %fwhm(i,n) = std(nVector);
                        %x0(i,n) = mean(nVector);
                        
                    catch
                        fwhm(i,n) = 1;
                        amplitude(i,n) = 0;
                        'try error 1'
                        
                    end
                    
                end
                
            end
            
        end
        
        try
            
            [a,b] = min(fwhm(observe(:,n) >= 0.9*max(observe(:,n)),n));
            
            index = find(fwhm(:,n) <= 1.25*a & fwhm(:,n) > 0);
            
            [~,minFwhmIndex] = min(fwhm(index,n));
            
            smallestIndex = index(minFwhmIndex);
            
            if a ~= 0
                best2(:,n) = vertcat(allNeurons.eletParameters(:,smallestIndex),x0(smallestIndex,n),fwhm(smallestIndex,n));
                
            else
                best2(:,n) = [0 0 0 0];
                'a == 0 error'
                
            end
            
        catch
            best2(:,n) = [0 0 0 0];
            'try error 2'
            
        end
        
    n
    
end