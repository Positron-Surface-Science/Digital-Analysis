i=0;

vectors = [];
responses = [];

vectorsC = cell(117,1);
vectorsCi = cell(95,1);
vectorsCt = cell(14,1);

while i <= 38
    
    ans1 = net(pulsesFront(:,i*5000+1:i*5000+5000));
    ans2 = net(pulsesRight(:,i*5000+1:i*5000+5000));
    ans3 = net(pulsesBack(:,i*5000+1:i*5000+5000));
    
    totalVector1 = sum(ans1,2);
    totalVector2 = sum(ans2,2);
    totalVector3 = sum(ans3,2);
    
    vectors = horzcat(vectors,totalVector1,totalVector2,totalVector3);
    
    vectorsC{i*3+1} = totalVector1;
    vectorsC{i*3+2} = totalVector2;
    vectorsC{i*3+3} = totalVector3;
    
    responses = horzcat(responses,0);
    %responsesC{i*3+1} = zeros(1,3);
    %responsesC{i*3+1}(1) = 1;
    responsesC{i*3+1} = 'front';
    %responsesC{i*3+1}(:) = num2str(responsesC{i*3+1}(:));
    
    responses = horzcat(responses,1);
    %responsesC{i*3+2} = zeros(1,3);
    %responsesC{i*3+2}(2) = 1;
    responsesC{i*3+2} = 'right';
    
    responses = horzcat(responses,2);
    %responsesC{i*3+3} = zeros(1,3);
    %responsesC{i*3+3}(3) = 1;
    responsesC{i*3+3} = 'back';
    
    i = i + 1;
    
end

for i=1:95
    
    vectorsCi{i}(:,1) = vectorsC{i}(:,1);
    %responsesCi{i}(:) = responsesC{i}(:);
    
end

for i=96:117
    
    vectorsCt{i-95}(:) = vectorsC{i}(:);
    %responsesCt{i-95}(:) = responsesC{i}(:);
    
end

responsesCiC = categorical(responsesCi);
responsesCtC = categorical(responsesCtC);