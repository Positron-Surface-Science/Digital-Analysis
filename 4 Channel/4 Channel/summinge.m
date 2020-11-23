
summer1 = zeros(512,1);
summer2 = zeros(512,1);

for E1 = 1:512
    for E2 = 1:512
        
            summer1(E1) = summer1(E1) + eVse(E1,E2);
            summer2(E2) = summer1(E2) + eVse(E1,E2);
            
    end
end

figure
plot([summer1, summer2]);