
summer = zeros(512,1);
for Ge = 1813:2048
    for ToF = 1:512
        %for Ge = 2633:3000
            summer(ToF) = summer(ToF) + gammaVsToF(Ge,ToF);
        %end
    end
end

summer = summer;
%sum = smooth(sum,25,'lowess');
figure
plot(summer);