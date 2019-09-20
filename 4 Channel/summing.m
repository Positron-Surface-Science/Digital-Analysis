
summer = zeros(1500,1);
for Ge = 1:1500
    for ToF = 1:1000
        %for Ge = 2633:3000
            summer(Ge) = summer(Ge) + gammaVsToF(Ge,ToF);
        %end
    end
end

summer = summer;
%sum = smooth(sum,25,'lowess');
figure
plot(summer);