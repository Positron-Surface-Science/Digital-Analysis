pBegin = 350
pEnd = 2000

r1Begin = 1
r1End = 2500

r2Begin = 2501
r2End = 5000

ps1 = horzcat(p1(pBegin:pEnd,r1Begin:r1End),p1CopyR(pBegin:pEnd,r1Begin:r1End),p1Copy1R(pBegin:pEnd,r1Begin:r1End),...p1Copy2R(pBegin:pEnd,r1Begin:r1End),p1Copy3R(pBegin:pEnd,r1Begin:r1End),p1Copy4R(pBegin:pEnd,r1Begin:r1End), ...
    p2(pBegin:pEnd,r1Begin:r1End),p2CopyR(pBegin:pEnd,r1Begin:r1End),p2Copy1R(pBegin:pEnd,r1Begin:r1End),...p2Copy2R(pBegin:pEnd,r1Begin:r1End),p1Copy3R(pBegin:pEnd,r1Begin:r1End),p2Copy4R(pBegin:pEnd,r1Begin:r1End), ...
    p3(pBegin:pEnd,r1Begin:r1End),p3CopyR(pBegin:pEnd,r1Begin:r1End),p3Copy1R(pBegin:pEnd,r1Begin:r1End),...p3Copy2R(pBegin:pEnd,r1Begin:r1End),p1Copy3R(pBegin:pEnd,r1Begin:r1End),p3Copy4R(pBegin:pEnd,r1Begin:r1End), ...
    p4(pBegin:pEnd,r1Begin:r1End),p4CopyR(pBegin:pEnd,r1Begin:r1End),p4Copy1R(pBegin:pEnd,r1Begin:r1End),...p4Copy2R(pBegin:pEnd,r1Begin:r1End),p1Copy3R(pBegin:pEnd,r1Begin:r1End),p4Copy4R(pBegin:pEnd,r1Begin:r1End), ...
    p5(pBegin:pEnd,r1Begin:r1End),p5CopyR(pBegin:pEnd,r1Begin:r1End),p5Copy1R(pBegin:pEnd,r1Begin:r1End));%,p5Copy2R(pBegin:pEnd,r1Begin:r1End),p1Copy3R(pBegin:pEnd,r1Begin:r1End),p5Copy4R(pBegin:pEnd,r1Begin:r1End));%, ...
    %pAs(pBegin:pEnd,1:1250));

ps2 = horzcat(p1(pBegin:pEnd,r2Begin:r2End),p1CopyR(pBegin:pEnd,r2Begin:r2End),p1Copy1R(pBegin:pEnd,r2Begin:r2End),...p1Copy2R(pBegin:pEnd,r2Begin:r2End),p1Copy3R(pBegin:pEnd,r2Begin:r2End),p1Copy4R(pBegin:pEnd,r2Begin:r2End), ...
    p2(pBegin:pEnd,r2Begin:r2End),p2CopyR(pBegin:pEnd,r2Begin:r2End),p2Copy1R(pBegin:pEnd,r2Begin:r2End),...p2Copy2R(pBegin:pEnd,r2Begin:r2End),p2Copy3R(pBegin:pEnd,r2Begin:r2End),p2Copy4R(pBegin:pEnd,r2Begin:r2End), ...
    p3(pBegin:pEnd,r2Begin:r2End),p3CopyR(pBegin:pEnd,r2Begin:r2End),p3Copy1R(pBegin:pEnd,r2Begin:r2End),...p3Copy2R(pBegin:pEnd,r2Begin:r2End),p3Copy3R(pBegin:pEnd,r2Begin:r2End),p3Copy4R(pBegin:pEnd,r2Begin:r2End), ...
    p4(pBegin:pEnd,r2Begin:r2End),p4CopyR(pBegin:pEnd,r2Begin:r2End),p4Copy1R(pBegin:pEnd,r2Begin:r2End),...p4Copy2R(pBegin:pEnd,r2Begin:r2End),p4Copy3R(pBegin:pEnd,r2Begin:r2End),p4Copy4R(pBegin:pEnd,r2Begin:r2End), ...
    p5(pBegin:pEnd,r2Begin:r2End),p5CopyR(pBegin:pEnd,r2Begin:r2End),p5Copy1R(pBegin:pEnd,r2Begin:r2End));%,p5Copy2R(pBegin:pEnd,r2Begin:r2End),p5Copy3R(pBegin:pEnd,r2Begin:r2End),p5Copy4R(pBegin:pEnd,r2Begin:r2End));%, ...
    %pAs(pBegin:pEnd,1251:2500));

ps1Baseline = 0;%mean(ps1(1:350,:));
ps2Baseline = 0;%mean(ps2(1:350,:));

ps1C = ps1 - ps1Baseline;
ps2C = ps2 - ps2Baseline;

    
ps1Mean = mean(mean(ps1C,1));
ps1Std = std(std(ps1C,1));

cData1 = num2cell((ps1C - ps1Mean)/ps1Std, 1);
cResponse1 = num2cell((ps2C - ps1Mean)/ps1Std, 1);
