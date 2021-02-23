%{ 
 Best so far:
 numHiddenUnits = 250;
 dropoutLayer(0.95, 'Name', 'dpo1')
    'InitialLearnRate',0.001, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.75, ...
    'LearnRateDropPeriod',50, ...
    'L2Regularization',0.02, ...

%}

numFeatures = 64;
numHiddenUnits = 10;
numResponses = 3;

layers = [ ...
    sequenceInputLayer(numFeatures, 'Name', 'input')
    
    %sequenceFoldingLayer('Name', 'fold')
    
    %convolution2dLayer([3 1], 8, 'Name', 'cnl1')
    %maxPooling2dLayer([4 1], 'Name', 'mp1')
    %fullyConnectedLayer(numHiddenUnits, 'Name', 'fc1')
    
    %leakyReluLayer('Name', 'lrelu1')
    bilstmLayer(numHiddenUnits+150, 'Name', 'lstm0', 'OutputMode','last')
    
    
    %sequenceUnfoldingLayer('Name', 'unfold')
    
    flattenLayer('Name', 'fll1')
    
    %lstmLayer(numHiddenUnits+1, 'Name', 'lstm1', 'OutputMode','last')
    
    %fullyConnectedLayer(numHiddenUnits, 'Name', 'fc1')
    reluLayer('Name', 'relu1')
    %dropoutLayer(0.1, 'Name', 'dpo1')
    
    fullyConnectedLayer(numResponses, 'Name', 'fc2')
    reluLayer('Name', 'relu2')
    
    
    softmaxLayer('Name', 'sml1')
    
    classificationLayer('Name', 'cl1')];

    %regressionLayer('Name', 'rl1')];

maxEpochs = 5000;
miniBatchSize = 64;

lgraph = layerGraph(layers);
%lgraph = addLayers(lgraph,sLayer1);
%lgraph = addLayers(lgraph,sLayer2);
%lgraph = connectLayers(lgraph,'relu','resizingConv1');
%lgraph = connectLayers(lgraph,'lrelu2','resizingConv2');
%lgraph = connectLayers(lgraph,'fold/miniBatchSize','unfold/miniBatchSize');
%lgraph = connectLayers(lgraph,'skippingBn1','add1/in2');
%lgraph = connectLayers(lgraph,'skippingBn2','add2/in2');


options = trainingOptions('adam', ...
    'ExecutionEnvironment','cpu', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',Inf, ...
    'ValidationFrequency',150, ...
    'InitialLearnRate',0.001, ...
    'ValidationData',{vectorsCt,responsesCt}, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.75, ...
    'LearnRateDropPeriod',550, ...
    'L2Regularization',0.0010, ...
    'Shuffle','every-epoch', ...
    'Verbose',true, ...
    'Plots','training-progress');

    %
net4 = trainNetwork(vectorsCi, responsesCi, lgraph, options);

%ans = classify(net4,testData);

%acc = sum(ans == testResponse)./numel(testResponse)