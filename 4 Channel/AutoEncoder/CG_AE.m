function [y_train, W, Wo] = CG_AE(x,T,Nh,Nit,W)

Nv = size(x,1);
X = [x ones(Nv,1)];
N = size(x,2);
M = N;
%-------A: NET CONTROL OF W INITIALIZATION-------

%Random Generation of weight matrix(Nh x N+1) with STD=1,Mean=0

W= padarray(W,[1,0],0,'post');

%-------B: OWO FOR OUTPUT WEIGHT (Wo) INITIALIZATION-------

OpNew = X*W';
OpNew(:,end) = 1 ;


R=( OpNew' * OpNew)/Nv;
C=( OpNew' * T)/Nv;

Wo = zeros(Nh+1, M);
for i = 1:M
    [Wo(:, i)] = CG_OWO(R, C(:,i));
end

Wo=Wo';


%-------C: GRADIENT CALCULATIONS--------
Train_Error = zeros(1,Nit);

Egp = 1 ;

Pa = zeros(Nh+1,N+1);
Po = zeros(M,Nh+1);

Nw = (Nh+1)*M;
for iter=1:Nit
    
     if floor(iter/Nw) * Nw == iter
        Pa = zeros(Nh+1,N+1);
        Po = zeros(M,Nh+1);
        Egp = 1;
    end
    
    
    OpNew = X*W';
    
    OpNew(:,end) =1 ;
    
    y = OpNew * Wo';
    
    ty= T-y;
    
    ga = Wo' * ( ty' * X ) ;
    ga=ga/Nv;
    
    if iter == 1
        go = zeros(M,Nh+1);
    else
        go = ty' * OpNew;
        go=go/Nv;
    end 
    
    %% Direction Vector

    Egi = sum(sum(ga .*ga));
    
    Ego = sum(sum(go .^2));
    
    EG = Egi + Ego;

    dir = EG / Egp;
    Egp = EG;

    Pa = ga + dir*Pa;
    
    Po = go + dir*Po;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% ------- D: OPTIMAL LEARNING FACTOR CALCULATION-------
    gx = Pa*X';
    D1 = 0;
    D2 = 0;
    dy = Wo*gx + Po* OpNew' ;
    
    for p = 1:Nv
        for i=1:M
            D1=D1+(T(p,i)-y(p,i))*dy(i,p);
            D2=D2+ dy(i,p)*dy(i,p);
        end
    end
    z = D1/D2;
    
    W = W + z*Pa ;
    
    Wo = Wo + z*Po ;
    
    
    % training ends here
    %---------------------------------------------%
    
    Opp = X*W' ;
    Opp(:,end) = 1;
    
    y_train=Opp*Wo';
    
    Ev=0;
    
    for i=1:N
        for p=1:Nv
            Ev =Ev + (x(p,i)-y_train(p,i))^2;
        end
    end
    
    Er_train =  Ev/Nv;
    Train_Error(1,iter) =Er_train;
    
    fprintf('\nIteration no: %d', iter)
    fprintf('\nTraining Error = %.8f \n', Er_train)

end
%--------- end of iterations -------------%


[MinTrainError,minIndexTrain] = min(Train_Error);

fprintf("\nMinimum training error = %f  at It = %d\n",MinTrainError, minIndexTrain)

 figure
 plot(Train_Error)
legend('Training Error')
end





