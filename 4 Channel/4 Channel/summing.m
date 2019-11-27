
summer = zeros(1024,1);
for Ge = 172:246
    for ToF = 1:1024
        %for Ge = 2633:3000
            summer(ToF) = summer(ToF) + electronVsElectron(Ge,ToF);
        %end
    end
end

summer = summer;
%sum = smooth(sum,25,'lowess');
figure
plot(summer);