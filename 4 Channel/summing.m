
summer = zeros(4096,1);
for Ge = 1:4096
    for ToF = 411:473
        %for Ge = 2633:3000
            summer(Ge) = summer(Ge) + gammaVsToF(Ge,ToF);
        %end
    end
end

summer = summer;
%sum = smooth(sum,25,'lowess');
figure
plot(summer);