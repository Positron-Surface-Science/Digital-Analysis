xToF = linspace(-0.25E-6,0.25E-6,1024)';


    observe = [];
    fwhm = [];
    index = [];
    best2 = [];
    

for n=1:400
    
    %if exact(n) == 1
        
        for i=1:260
            
            if allNeurons.eletParameters(1,i) ~= 0
                
                allNeurons.neurons{n}(allNeurons.neurons{n}(:,i) == 0,i) = NaN;
                
                hToF = histcounts(allNeurons.neurons{n}(:,i),'NumBins',1024,'BinLimits',[-0.25E-6 0.25E-6])';
                hToF = hToF/max(hToF);
                
                observe(i,n) = numel(allNeurons.neurons{n}(isnan(allNeurons.neurons{n}(:,i)) == 0,i));
                
                if max(hToF) < sum(hToF) && observe(i,n) >= 0.5*max(observe(:,n))
                    
                    try
                        
                        [~,second] = max(hToF);
                        
                        gaussFit = fit(xToF(second-125:second+125),hToF(second-125:second+125),'gauss1');
                        
                        cValues = coeffvalues(gaussFit);
                        
                        
                        fwhm(i,n) = (2*sqrt(log(2))*cValues(3))/sqrt(2);
                        amplitude(i,n) = cValues(1);
                        x0(i,n) = cValues(2);
                        
                        
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
                        
                    end
                    
                end
                
            end
            
        end
        
        try
            
            [a,b] = min(fwhm(observe(:,n) >= 0.5*max(observe(:,n)),n));
            
            index = find(fwhm(:,n) <= 1.05*a & fwhm(:,n) > 0);
            
            [~,minFwhmIndex] = min(fwhm(index,n));
            
            smallestIndex = index(minFwhmIndex);
            
            if a ~= 0
                best2(:,n) = vertcat(allNeurons.eletParameters(:,smallestIndex),x0(smallestIndex,n),fwhm(smallestIndex,n));
                'sdgsdaf'
                
            else
                best2(:,n) = [0 0 0 0];
                
            end
            
        catch
            best2(:,n) = [0 0 0 0];
            
        end
        
    n
    
end