function [Wo] = CG_OWO(R, Ci)

    Nu = size(R, 1);
    Wo = zeros(Nu, 1);

		p = zeros(Nu, 1);
		
		XD = 1;
		% Loop over iterations
        for iter = 1:Nu
            % Calculate gradient
            g = -2 * (Ci - R * Wo);
            
            % Calculate enrgy in the gradient and B1
            XN = sum(g .* g);
            B1 = XN / XD;
            XD = XN;
    
            % Find direction vector
            p = -g + B1 * p;
            
            % Calculate B2
            Num = -sum(p .* g)/2;
            Den = p' * R * p;
            B2 = Num / Den;
            
           
            Wo = Wo + B2 * p;
        end
end