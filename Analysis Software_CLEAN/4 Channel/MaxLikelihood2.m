format long
digits(15)
number = 174538;

pulses = horzcat(pIn,pOut);

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

for it=1:10
    
    % Solving for a
    
    for i=1:number
        
        sum1 = 0;
        sum2 = 0;
        
        for n=1:numel(pulses(:,i)-1)
            
            sum1 = sum1 + pulses(n,i)*s(n) + pulses(n,i)*b(i);
            sum2 = sum2 + b(i)^2 + s(n)^2 + 2*b(i)*s(n);
            
        end
        
        a(i) = sum1/sum2;
        
    end
    
    % Solving for b
    
    for i=1:number
        
        sum1 = 0;
        sum2 = 0;
        
        for n=1:numel(pulses(:,i)-1)
            
            sum1 = sum1 + pulses(n,i)*a(i) - a(i)^2*s(n);
            
        end
        
        sum2 = a(i)^2;
        
        b(i) = sum1/sum2;
        
    end
    
    it
    
    % Solving for s
    
    for n=1:numel(pulses(:,i)-1)
        
        sum1 = 0;
        sum2 = 0;
        
        for i=1:number
            
            sum1 = sum1 + a(i)*pulses(n,i) -  b(i)*a(i)^2;
            sum2 = sum2 + a(i)^2;
            
        end
        
        s(n) = sum1/sum2;
        
    end
    
end

sum2 = 0;

% Calculating x(i)

x = zeros(size(pulses));

for i=1:number
    
    x(:,i) = a(i)*s - b(i);
    
    err(i) = immse(x(:,i),pulses(:,i));
    
end


