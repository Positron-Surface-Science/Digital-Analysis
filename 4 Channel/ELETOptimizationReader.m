ELETMatrix = zeros(length(exact),3);
noFolder = false;

try
    folderSelections = 'C:\Saved Data Analysis - Lecroy\Ge\2019\September\HPGe_NaI_Coincidence_20kV_Zeolite-AsReceived\ANNDataCollection\Neurons'%uigetdir;
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
                nif = fileList(t).name;
                e = strfind(nif,'_');
                p = strfind(nif,'%');
                k = strfind(nif,'FWHM');
                
                if isempty(k) == 0
                    numbers = nif(e(5)+1:k-1);
                    values(t) = str2double(numbers);
                    
                    if t > 1 && values(t) < values(t-1) && ...
                            values(t) < newValue
                        elet1N = nif(e(3)+1:p(1)-1);
                        elet2N = nif(e(4)+1:p(2)-1);
                        
                        ELETVector(neuron,1) = str2double(elet1N);
                        ELETVector(neuron,2) = str2double(elet2N);
                        newValue = values(t);
                        
                        goodFile(n) = t;
                        
                    end
                    
                end
            end
            
            if goodFile(n) > 0
                goodMatrix = importdata([folderSelections,'\',folderList(n).name,'\',fileList(goodFile(n)).name]);
                
                gaussFit = fit(x,goodMatrix,'gauss1');
                
                cValues = coeffvalues(gaussFit);
                
                ELETVector(neuron,3) = cValues(2);
            end
    end
    
end




                    