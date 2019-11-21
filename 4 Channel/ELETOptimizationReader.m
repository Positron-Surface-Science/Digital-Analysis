ELETMatrix = zeros(length(exact),4);
noFolder = false;
FWHMVector = zeros(length(exact),1);

try
    folderSelections = 'C:\Saved Data Analysis - Lecroy\Ge\2019\September\HPGe_BaF2_Coincidence_1000V_ZeolitePowder-Baked\ANN Analysis On\9x25Network\Neurons'%uigetdir;
catch
    noFolder = true;
end

if noFolder == false && isnumeric(folderSelections) == false && exist(folderSelections,'dir') ~= 0
    folderList = dir(folderSelections);
    goodFile = zeros(numel(folderList),1);
    
    for n=3:numel(folderList)
        fileList = dir([folderSelections,'\',folderList(n).name,'\*.txt']);
        numFiles = length(fileList);
        values = zeros(numFiles,1);
        folderName = folderList(n).name;
        neuronFind = strfind(folderName,' ');
        neuronS = folderName(neuronFind+1:length(folderName));
        neuron = str2double(neuronS);
        newValue = 1;
        
        for t=1:numFiles
            %try
                nif = fileList(t).name;
                e = strfind(nif,'_');
                p = strfind(nif,'%');
                k = strfind(nif,'FWHM');
                
                if isempty(k) == 0
                    n
                    goodMatrix = importdata([folderSelections,'\',folderName,'\',nif]);
                    
                    gaussFit = fit(xdataToF,goodMatrix,'gauss1');
                    
                    cValues = coeffvalues(gaussFit);
                    
                    numbers = (2*sqrt(log(2))*cValues(3))/sqrt(2);%nif(e(5)+1:k-1);
                    values(t) = numbers;
                    
                    if t > 1 && values(t) < values(t-1) && ...
                            values(t) < newValue && ...
                            values(t) > 0
                        elet1N = nif(e(3)+1:p(1)-1);
                        elet2N = nif(e(4)+1:p(2)-1);
                        
                        ELETMatrix(neuron,1) = str2double(elet1N);
                        ELETMatrix(neuron,2) = str2double(elet2N);
                        ELETMatrix(neuron,3) = cValues(2);
                        ELETMatrix(neuron,4) = values(t);
                        
                        newValue = values(t);
                        
                        goodFile(n) = t;
                        
                    end
                    
                end
           % catch
            %end
        end
        
        %{
        if goodFile(n) > 0
            goodMatrix = importdata([folderSelections,'\',folderList(n).name,'\',fileList(goodFile(n)).name]);
            
            gaussFit = fit(x,goodMatrix,'gauss1');
            
            cValues = coeffvalues(gaussFit);
            
            ELETMatrix(neuron,3) = cValues(2);
            
        end
        %}
    end
    
end




