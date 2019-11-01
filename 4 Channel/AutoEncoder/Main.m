%clc
%{
close all
%clear

format = input('Enter 1 for .tra format files, 2 for .csv files: ');

% %% For files ending with .tra : 
if format == 1 
    
    training_file = input('Enter file name: ', 's');
    N= input('Enter the number of Inputs: ');
    M = input('Enter number of Outputs: ');
    
    [x, t, Nv] = read_Approximation(training_file, N, M);

elseif format == 2 

    %% For .csv files  :
    training_file = input('Enter file name: ', 's');
%}
    N= input('Enter the number of Inputs: ');
    M = input('Enter number of Outputs: ');
  %{  
    file = load(training_file);
    
    x = file( : , 1: N);
end
%}
%% 
T = x; % setting target equal to inputs.

Nh= input('Enter the number of Hidden nodes: ');
Nit = input('Enter number of iterations: ');

rng(1)
W= normrnd(0,1,Nh,N+1) ;

[y_train, W1, Wo1] = CG_AE(x,T,Nh,Nit,W);

%% Testing code  (Uncomment this section if you  want to check testing Pe)
% 
% %% For files ending with .tra : 
% if format == 1 
%     
%     testing_file = input('\nEnter test file name: ', 's');
%    
%     [x_test, t_test, Nv_test] = read_Approximation(testing_file, N, M);
% 
% elseif format == 2 
% 
%     %% For .csv files  :
%     testing_file = input('Enter test file name: ', 's');
%     file_test = load(testing_file);
%     
%     x_test = file_test( : , 1: 8);
%     Nv_test = size(x_test,1);
% 
% end
% 
% %%
% X_test = [x_test ones(Nv_test,1)];
% O_test = X_test * W1';
% O_test(:,end) = 1;
% 
% y_test = O_test * Wo1';
% 
% Ev_test=0;
% 
% for i=1:N
%     for p=1:Nv_test
%         Ev_test =Ev_test + (x_test(p,i)-y_test(p,i))^2;
%     end
% end
% 
% Er_test =  Ev_test/Nv_test;
% 
% fprintf('\nTestError = %.8f \n',Er_test)

  