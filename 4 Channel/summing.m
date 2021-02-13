
summer = zeros(2048,1);
for Ge = 1:2048
    for ToF = 81:134
        %for Ge = 2633:3000
            summer(Ge) = summer(Ge) + gammaVsToF(Ge,ToF);
        %end
    end
end

summer = summer;
%sum = smooth(sum,25,'lowess');
figure
plot(summer);