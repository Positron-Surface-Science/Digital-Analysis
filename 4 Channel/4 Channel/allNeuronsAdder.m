%{

function allNeurons = allNeuronsAdder(allNeurons1, allNeurons2, allNeurons3, ...
    allNeurons4, allNeurons5, allNeurons6, allNeurons7, allNeurons8, allNeurons9)

for c = 1:9
    
    eval(['aElet = ', 'allNeurons', c, '.eletParameters(:,1:50);']);
    
    newC = 50*c - 49;
    
    allNeurons.eletParameters(:,newC:newC + 50) = aElet;
    
    for y = 1:1024
        
       eval(['aNeurons = ', 'allNeurons', c, '.neurons(', y, ')(:,1:50);']);
        
       allNeurons.neurons{y}(:,newC:newC + 50) = aNeurons; 
        
    end
    
end

end

%}


allNeurons.eletParameters(:,1:50) = allNeurons1.eletParameters(:,1:50);

for y=1:1024
allNeurons.neurons{y}(:,1:50) = allNeurons1.neurons{y}(:,1:50);
end

allNeurons.eletParameters(:,51:100) = allNeurons2.eletParameters(:,1:50);

for y=1:1024
allNeurons.neurons{y}(:,51:100) = allNeurons2.neurons{y}(:,1:50);
end

allNeurons.eletParameters(:,101:150) = allNeurons3.eletParameters(:,1:50);

for y=1:1024
allNeurons.neurons{y}(:,101:150) = allNeurons3.neurons{y}(:,1:50);
end

allNeurons.eletParameters(:,151:200) = allNeurons4.eletParameters(:,1:50);

for y=1:1024
allNeurons.neurons{y}(:,151:200) = allNeurons4.neurons{y}(:,1:50);
end

allNeurons.eletParameters(:,201:250) = allNeurons5.eletParameters(:,1:50);

for y=1:1024
allNeurons.neurons{y}(:,201:250) = allNeurons5.neurons{y}(:,1:50);
end

allNeurons.eletParameters(:,251:300) = allNeurons6.eletParameters(:,1:50);

for y=1:1024
allNeurons.neurons{y}(:,251:300) = allNeurons6.neurons{y}(:,1:50);
end

allNeurons.eletParameters(:,301:350) = allNeurons7.eletParameters(:,1:50);

for y=1:1024
allNeurons.neurons{y}(:,301:350) = allNeurons7.neurons{y}(:,1:50);
end

allNeurons.eletParameters(:,351:400) = allNeurons8.eletParameters(:,1:50);

for y=1:1024
allNeurons.neurons{y}(:,351:400) = allNeurons8.neurons{y}(:,1:50);
end

allNeurons.eletParameters(:,401:450) = allNeurons9.eletParameters(:,1:50);

for y=1:1024
allNeurons.neurons{y}(:,401:450) = allNeurons9.neurons{y}(:,1:50);
end
