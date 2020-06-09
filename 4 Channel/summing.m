
summer = zeros(2048,1);
for Ge = 373:674
    for ToF = 1:2048
        %for Ge = 2633:3000
            summer(Ge) = summer(Ge) + electronVsElectron(Ge,ToF);
        %end
    end
end

summer = summer;
%sum = smooth(sum,25,'lowess');
figure
plot(summer);