
summer = zeros(1024,1);
for Ge = 1:1024
    for ToF = 195:1035
        %for Ge = 2633:3000
            summer(Ge) = summer(Ge) + gammaVsToF(Ge,ToF);
        %end
    end
end

summer = summer;
%sum = smooth(sum,25,'lowess');
figure
plot(summer);