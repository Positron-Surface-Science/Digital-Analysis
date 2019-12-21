% Solve a Clustering Problem with a Self-Organizing Map
% Script generated by Neural Clustering app
% Created 07-Oct-2019 16:47:58
%
% This script assumes these variables are defined:
%
%   p5s - input data.

x = nT(:,1:75000);

% Create a Self-Organizing Map
dimension1 = 20;
dimension2 = 20;
net = selforgmap([dimension1 dimension2],150,4,'hextop','linkdist');

% Choose Plot Functions
% For a list of all plot functions type: help nnplot
%net.plotFcns = {'plotsomtop','plotsomnc','plotsomnd', ...
    %'plotsomplanes', 'plotsomhits', 'plotsompos'};

net.trainParam.epochs = 2000;

% Train the Network
[net,tr] = train(net,x);

% Test the Network
y = net(x);

% View the Network
view(net)

% Plots
% Uncomment these lines to enable various plots.
%figure, plotsomtop(net)
%figure, plotsomnc(net)
%figure, plotsomnd(net)
%figure, plotsomplanes(net)
%figure, plotsomhits(net,x)
%figure, plotsompos(net,x)

% Deployment
% Change the (false) values to (true) to enable the following code blocks.
% See the help for each generation function for more information.
if (false)
    % Generate MATLAB function for neural network for application
    % deployment in MATLAB scripts or with MATLAB Compiler and Builder
    % tools, or simply to examine the calculations your trained neural
    % network performs.
    genFunction(net,'myNeuralNetworkFunction');
    y = myNeuralNetworkFunction(x);
end
if (false)
    % Generate a matrix-only MATLAB function for neural network code
    % generation with MATLAB Coder tools.
    genFunction(net,'myNeuralNetworkFunction','MatrixOnly','yes');
    y = myNeuralNetworkFunction(x);
end
if (false)
    % Generate a Simulink diagram for simulation or deployment with.
    % Simulink Coder tools.
    gensim(net);
end
