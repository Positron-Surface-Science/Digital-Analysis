
summer = zeros(512,1);

for Ge = 143:256
    for ToF1 = 1:512
        for ToF2 = 1:512
            %for Ge = 2633:3000
            summer(ToF2) = summer(ToF2) + gammaVsToF(Ge,ToF1,ToF2);
            %end
        end
    end
end

summer = summer;
%sum = smooth(sum,25,'lowess');
figure
plot(summer);