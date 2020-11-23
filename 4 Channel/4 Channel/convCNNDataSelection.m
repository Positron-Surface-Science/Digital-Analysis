

p1S = single(p1s);
p2S = single(p2s);
p3S = single(p3s);
p4S = single(p4s);
p5S = single(p5s);

%{
b1 = mean(p1s(200:400,:));
p1S = p1S - b1;
b2 = mean(p2s(200:400,:));
p2S = p2S - b2;
b3 = mean(p3s(200:400,:));
p3S = p3S - b3;
b4 = mean(p4s(200:400,:));
p4S = p4S - b4;
b5 = mean(p5s(200:400,:));
p5S = p5S - b5;
%}
pBegin = 501
pEnd = 4000

r1Begin = 1
r1End = 7000

r2Begin = 1
r2End = 3500

E1 = r1End+11000;
E2 = r1End;
E3 = r1End+10000;
E4 = r1End+39000;
E5 = r1End-1000;

%{
cols = size(p1S, 2);
P = randperm(cols);
p1CopyR = p1S(:, P);

cols = size(p1S, 2);
P = randperm(cols);
p1Copy1R = p1S(:, P);

cols = size(p2S, 2);
P = randperm(cols);
p2CopyR = p2S(:, P);

cols = size(p2S, 2);
P = randperm(cols);
p2Copy1R = p2S(:, P);

cols = size(p3S, 2);
P = randperm(cols);
p3CopyR = p3S(:, P);

cols = size(p3S, 2);
P = randperm(cols);
p3Copy1R = p3S(:, P);

cols = size(p4S, 2);
P = randperm(cols);
p4CopyR = p4S(:, P);

cols = size(p4S, 2);
P = randperm(cols);
p4Copy1R = p4S(:, P);

cols = size(p5S, 2);
P = randperm(cols);
p5CopyR = p5S(:, P);

cols = size(p5S, 2);
P = randperm(cols);
p5Copy1R = p5S(:, P);

%}

ps1 = horzcat(p1S(pBegin:pEnd,r1Begin:r1End+11000),...p1CopyR(pBegin:pEnd,r1Begin:r1End),p1Copy1R(pBegin:pEnd,r1Begin:r1End),...p1Copy2R(pBegin:pEnd,r1Begin:r1End),p1Copy3R(pBegin:pEnd,r1Begin:r1End),p1Copy4R(pBegin:pEnd,r1Begin:r1End), ...
    p2S(pBegin:pEnd,r1Begin:r1End),...p2CopyR(pBegin:pEnd,r1Begin:r1End),p2Copy1R(pBegin:pEnd,r1Begin:r1End),...p2Copy2R(pBegin:pEnd,r1Begin:r1End),p1Copy3R(pBegin:pEnd,r1Begin:r1End),p2Copy4R(pBegin:pEnd,r1Begin:r1End), ...
    p3S(pBegin:pEnd,r1Begin:r1End+10000),...p3CopyR(pBegin:pEnd,r1Begin:r1End),p3Copy1R(pBegin:pEnd,r1Begin:r1End),...p3Copy2R(pBegin:pEnd,r1Begin:r1End),p1Copy3R(pBegin:pEnd,r1Begin:r1End),p3Copy4R(pBegin:pEnd,r1Begin:r1End), ...
    p4S(pBegin:pEnd,r1Begin:r1End+39000),...p4CopyR(pBegin:pEnd,r1Begin:r1End),p4Copy1R(pBegin:pEnd,r1Begin:r1End),...p4Copy2R(pBegin:pEnd,r1Begin:r1End),p1Copy3R(pBegin:pEnd,r1Begin:r1End),p4Copy4R(pBegin:pEnd,r1Begin:r1End), ...
    p5S(pBegin:pEnd,r1Begin:r1End-1000));%,p5CopyR(pBegin:pEnd,r1Begin:r1End),p5Copy1R(pBegin:pEnd,r1Begin:r1End));%,p5Copy2R(pBegin:pEnd,r1Begin:r1End),p1Copy3R(pBegin:pEnd,r1Begin:r1End),p5Copy4R(pBegin:pEnd,r1Begin:r1End));%, ...
    %pAs(pBegin:pEnd,1:1250));

    
    ps2 = zeros([1 37500]);
    
    ps2(r1Begin:E1) = 80.9979;
    ps2(E1+1:E1+E2) = 276.3992;
    ps2(E1+E2+1:E1+E2+E3) = 302.8512;
    ps2(E1+E2+E3+1:E1+E2+E3+E4) = 356.0134;
    ps2(E1+E2+E3+E4+1:E1+E2+E3+E4+E5) = 383.8491;
    
    %{
ps2 = horzcat(p1S(pBegin:pEnd,r2Begin:r2End+11000),...p1CopyR(pBegin:pEnd,r2Begin:r2End),p1Copy1R(pBegin:pEnd,r2Begin:r2End),...p1Copy2R(pBegin:pEnd,r2Begin:r2End),p1Copy3R(pBegin:pEnd,r2Begin:r2End),p1Copy4R(pBegin:pEnd,r2Begin:r2End), ...
    p2S(pBegin:pEnd,r2Begin:r2End),...p2CopyR(pBegin:pEnd,r2Begin:r2End),p2Copy1R(pBegin:pEnd,r2Begin:r2End),...p2Copy2R(pBegin:pEnd,r2Begin:r2End),p2Copy3R(pBegin:pEnd,r2Begin:r2End),p2Copy4R(pBegin:pEnd,r2Begin:r2End), ...
    p3S(pBegin:pEnd,r2Begin:r2End+10000),...p3CopyR(pBegin:pEnd,r2Begin:r2End),p3Copy1R(pBegin:pEnd,r2Begin:r2End),...p3Copy2R(pBegin:pEnd,r2Begin:r2End),p3Copy3R(pBegin:pEnd,r2Begin:r2End),p3Copy4R(pBegin:pEnd,r2Begin:r2End), ...
    p4S(pBegin:pEnd,r2Begin:r2End+39000),...p4CopyR(pBegin:pEnd,r2Begin:r2End),p4Copy1R(pBegin:pEnd,r2Begin:r2End),...p4Copy2R(pBegin:pEnd,r2Begin:r2End),p4Copy3R(pBegin:pEnd,r2Begin:r2End),p4Copy4R(pBegin:pEnd,r2Begin:r2End), ...
    p5S(pBegin:pEnd,r2Begin:r2End-1000));%,p5CopyR(pBegin:pEnd,r2Begin:r2End),p5Copy1R(pBegin:pEnd,r2Begin:r2End));%,p5Copy2R(pBegin:pEnd,r2Begin:r2End),p5Copy3R(pBegin:pEnd,r2Begin:r2End),p5Copy4R(pBegin:pEnd,r2Begin:r2End));%, ...
    %pAs(pBegin:pEnd,1251:2500));
%}

for j=1:length(ps2(1,:))
    
    measurePulse1 = smooth(ps1,51,'moving');
    %measurePulse2 = smooth(ps2,51,'moving');
    
    
    [m1,ending1] = max(measurePulse1);
    %[m2,ending2] = max(measurePulse2);
    
    [~,begin1] = min(abs(measurePulse1-0.1*m1));
    %[~,begin2] = min(abs(measurePulse2-0.1*m2));
    
    riseTime1 = ending1 - begin1;
    %riseTime2 = ending2 - begin2;
    
    if riseTime1 > 100 %|| riseTime2 > 100
        
        ps1(:,j) = [];
        ps2(:,j) = [];
        
    end
    
end
    
    
%ps1Baseline = mean(ps1(1:200,:));
%ps2Baseline = mean(ps2(1:200,:));

%ps1C = ps1 - ps1Baseline;
%ps2C = ps2 - ps2Baseline;

ps1Mean = mean(mean(ps1,1));
ps1Std = std(std(ps1,1));

cData = num2cell((ps1 - ps1Mean)/ps1Std, 1);
cResponse = num2cell((ps2 - ps1Mean)/ps1Std, 1);
