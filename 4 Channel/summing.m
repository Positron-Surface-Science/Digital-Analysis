
summer = zeros(2048,1);
for Ge = 1:1024
    for ToF = 1:2048
        %for Ge = 2633:3000
            summer(ToF) = summer(ToF) + gammaVsToF(Ge,ToF);
        %end
    end
end

summer = summer;
%sum = smooth(sum,25,'lowess');
figure
plot(summer);