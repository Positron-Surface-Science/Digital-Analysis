numFeatures1 = 3500;
numFeatures2 = 1;
numHiddenUnits1 = 4500;
numHiddenUnits2 = 1500;
numHiddenUnits3 = 1000;
numResponses = 1;

layers = [ ...
    
    sequenceInputLayer([numFeatures1 numFeatures2 1], 'Name', 'input')
    %lstmLayer(numHiddenUnits2, 'Name', 'bll1')
    sequenceFoldingLayer('Name', 'fold')
    
    % First Convolution
    convolution2dLayer([10 1], 16, 'Name', 'conv1', 'NumChannels', 1, 'Stride', 2, 'Padding', 'same')
    batchNormalizationLayer('Name', 'bn1')
    %reluLayer('Name', 'relu1')
    maxPooling2dLayer([3 1], 'Name', 'mpling1', 'Stride', 2)

    % Second Convolution
    %convolution2dLayer([5 9], 24, 'Name', 'conv2', 'NumChannels', 16, 'Stride', 1, 'Padding', 'same')
    %batchNormalizationLayer('Name', 'bn2')
    %leakyReluLayer('Name', 'relu2')
    %maxPooling2dLayer([2 2], 'Name', 'mpling2', 'Stride', 2)
    
    % Third Convolution
    %convolution2dLayer([3 1], 128, 'Name', 'conv3', 'NumChannels', 64, 'Stride', 1)
    %batchNormalizationLayer('Name', 'bn3')
    %tanhLayer('Name', 'relu3')
    %maxPooling2dLayer([2 1], 'Name', 'mpling3', 'Stride', 1)
    
    % Flatten Layer
    %flattenLayer('Name', 'flatten1')
    
    %fullyConnectedLayer(numHiddenUnits3, 'Name', 'fc1')
    %dropoutLayer(0.5, 'Name', 'dpo1')
    %tanhLayer('Name', 'relu1')
    
    %fullyConnectedLayer(numHiddenUnits2, 'Name', 'fc2')
    %dropoutLayer(0.5, 'Name', 'dpo2')
    %reluLayer('Name', 'relu2')
    
    % Fully Connected Layer
    fullyConnectedLayer(numHiddenUnits1, 'Name', 'fc3')
    dropoutLayer(0.5, 'Name', 'dpo3')
    %tanhLayer('Name', 'relu3')
    
    %fullyConnectedLayer(numHiddenUnits2, 'Name', 'fc4')
    %dropoutLayer(0.5, 'Name', 'dpo4')
    %reluLayer('Name', 'relu4')
    
    %fullyConnectedLayer(numHiddenUnits3, 'Name', 'fc5')
    %dropoutLayer(0.5, 'Name', 'dpo5')
    %tanhLayer('Name', 'relu5')
    
    % First Upsampling
    %transposedConv2dLayer([2 3], 24, 'Name', 'tconv1', 'Stride', 1)
    %leakyReluLayer('Name', 'relu8')
    %additionLayer(2, 'Name', 'add1')
    
    %fullyConnectedLayer(numHiddenUnits3, 'Name', 'fc8')
    %tanhLayer('Name', 'relu8')
    
    %fullyConnectedLayer(numHiddenUnits2, 'Name', 'fc9')
    %reluLayer('Name', 'relu9')
    
    %fullyConnectedLayer(numHiddenUnits1, 'Name', 'fc5')
    %leakyReluLayer('Name', 'relu7')
    
    % Second Upsampling
    %transposedConv2dLayer([3 5], 16, 'Name', 'tconv2', 'Stride', 1)
    %batchNormalizationLayer('Name', 'bn4')
    %leakyReluLayer('Name', 'relu3')
    
    
    
    % Third Upsampling
    %maxUnpooling2dLayer('Name', 'munpling')
    %transposedConv2dLayer(5, 32, 'Name', 'tconv3', 'Stride', 1)
    %batchNormalizationLayer('Name', 'bn5')
    %)
    
    %additionLayer(2, 'Name', 'add1')
    
    sequenceUnfoldingLayer('Name', 'unfold')
    %leakyReluLayer('Name', 'relu5')
    % Second Flatten Layer
    flattenLayer('Name', 'flatten2')
    
    %bilstmLayer(numHiddenUnits2, 'Name', 'bll1')
    
    
    %leakyReluLayer('Name', 'relu7')
    
    fullyConnectedLayer(numResponses, 'Name', 'fc7')
    %leakyReluLayer('Name', 'tanh1')
    %dropoutLayer(0.1, 'Name', 'dpo2')
    
    regressionLayer('Name', 'rl')];

maxEpochs = 2000;
miniBatchSize = 1024;

%sLayer1 = skipLayer(32,2,8,'1');
%sLayer2 = skipLayer(256,4,128,'2');

lgraph = layerGraph(layers);
lgraph = connectLayers(lgraph, 'fold/miniBatchSize', 'unfold/miniBatchSize');


%lgraph = addLayers(lgraph,sLayer1);
%lgraph = addLayers(lgraph,sLayer2);
%lgraph = connectLayers(lgraph,'relu1','resizingConv1');
%lgraph = connectLayers(lgraph,'lrelu2','resizingConv2');
%lgraph = connectLayers(lgraph, 'relu4', 'relu5');
%lgraph = connectLayers(lgraph,'skippingBn1','add1/in2');
%lgraph = connectLayers(lgraph,'skippingBn2','add2/in2');
%}

options = trainingOptions('adam', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'InitialLearnRate',0.01, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.75, ...
    'LearnRateDropPeriod',45, ...
    'GradientThreshold',Inf, ...
    'L2Regularization',0.0001, ...
    'Shuffle','every-epoch', ...
    'Plots','training-progress', ...
    'Verbose',1, ...
    'ExecutionEnvironment', 'cpu');

%'L2Regularization',0.1, ...

net4 = trainNetwork(cData, cResponse, lgraph, options);


function layers = convolutionalUnit(numF, stride, numC, tag)
layers = [ ...
    convolution2dLayer(5,numF,'NumChannels',numC,'Padding','same','Stride',stride,'Name',[tag,'conv1'])
    batchNormalizationLayer('Name',[tag,'BN1'])
    tanhLayer('Name',[tag,'elu1'])
    
    convolution2dLayer(5,numF,'NumChannels',numC,'Padding','same','Stride',stride,'Name',[tag,'conv2'])
    batchNormalizationLayer('Name',[tag,'BN2'])];
    
end

function skippingLayer = skipLayer(numF, stride, numC, tag)
skippingLayer = [ ...
    
    convolution2dLayer([1 1], numF, 'Name', ['resizingConv' tag], 'NumChannels', numC, 'Stride', stride, 'Padding', 4)
    batchNormalizationLayer('Name', ['skippingBn', tag])];

end
