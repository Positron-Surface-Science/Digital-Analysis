format long
digits(15)
pulses = horzcat(pIn,pOut);
number = numel(pulses(1,:));
a = [];
s = [];
b = [];
err = 0;
 
%{
for i = 1:number
    rng;
    x = rand;
    
    if a1(i)*100 <= 2.5 && x >= 0.3
        pulses(:,i) = 0;
        
    end
    
end
%}
%{
for i=1:number
    
undone = log(pulses(:,i));
r2 = r*max(undone);

undone(550:numel(undone)) = undone(550:numel(undone)) - r2;

pulses2(:,i) = real(exp(undone));

end
%}
% Initializing for a and b

for i=1:number
    
    sum1 = 0;
    
    for n=1:numel(pulses(:,1)-1)
        
        sum1 = sum1 + pulses(n,i)^2;
        
    end
    
    a(i) = sqrt(sum1/numel(pulses(:,i)));
    
    b(i) = 0;
    
end

% Initializing for s

for n=1:numel(pulses(:,i)-1)
    
    sum1 = 0;
    sum2 = 0;
    
    for i=1:number
        
        sum1 = sum1 + a(i)*pulses(n,i) - b(i)*a(i);
        sum2 = sum2 + a(i)^2;
        
    end
    
    s(n) = sum1/sum2;
    
end

sum1 = 0;
sum2 = 0;

% Iterating

for it=1:1
    
    % Solving for a
    
    for i=1:number
        
        sum1 = 0;
        sum2 = 0;
        
        for n=1:numel(pulses(:,i)-1)
            
            sum1 = sum1 + pulses(n,i)*s(n) - b(i)*s(n);
            sum2 = sum2 + s(n)^2;
            
        end
        
        a(i) = sum1/sum2;
        
    end
    
    % Solving for s
    
    for n=1:numel(pulses(:,i)-1)
        
        sum1 = 0;
        sum2 = 0;
        
        for i=1:number
            
            sum1 = sum1 + pulses(n,i)*a(i) - b(i)*a(i);
            sum2 = sum2 + a(i)^2;
            
        end
        
        s(n) = sum1/sum2;
        
    end
    
    % Solving for b
    
    for i=1:number
        
        sum1 = 0;
        
        for n=1:numel(pulses(:,i)-1)
            
            sum1 = sum1 + pulses(n,i) - a(i)*s(n);
            
        end
        
        b(i) = sum1/numel(pulses(:,i));
        
    end
    
    it
    
end

sum2 = 0;

% Calculating x(i)

x = zeros(size(pulses));

for i=1:number
    
    x(:,i) = a(i)*s + b(i);
    
    err(i) = immse(x(:,i),pulses(:,i));
    
end

totalerr = sum(err)

histogram(err)

