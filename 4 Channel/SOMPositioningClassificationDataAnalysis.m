i=0;

vectors = [];
responses = [];

vectorsC = cell(117,1);
vectorsCi = cell(95,1);
vectorsCt = cell(14,1);

numberInSet = 5000;

valSet = floor(0.1*195000/numberInSet);

while i <= floor(195000/numberInSet) - 1
    
    ans1 = net(pulsesFront(:,i*numberInSet+1:i*numberInSet+numberInSet));
    ans2 = net(pulsesRight(:,i*numberInSet+1:i*numberInSet+numberInSet));
    ans3 = net(pulsesBack(:,i*numberInSet+1:i*numberInSet+numberInSet));
    
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

%responsesCiC = categorical(responsesCi);
%responsesCtC = categorical(responsesCtC);

for i=1:95
vectorsCi{i}(:,1) = vectorsC{i}(:,1);
%responsesCt{i-95}(:) = responsesC{i}(:);
end

vectorsCi=vectorsCi';

for i=96:117
vectorsCt{i-95}(:,1) = vectorsC{i}(:,1);
%responsesCt{i-95}(:) = responsesC{i}(:);
end
vectorsCt=vectorsCt';

ans = categorical(responsesC);

responsesCi = ans(1:95);
responsesCt = ans(96:117);
responsesCi = responsesCi';
responsesCt = responsesCt';