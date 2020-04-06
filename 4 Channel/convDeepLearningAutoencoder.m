numFeatures = length(ps1(:,1));
numHiddenUnits1 = 50;
numHiddenUnits2 = 100;
numResponses1 = 500;
numResponses2 = 250;
numResponses3 = length(ps2(:,1));

layers = [ ...
    
    sequenceInputLayer([numFeatures 1 1], 'Name', 'input')
    sequenceFoldingLayer('Name', 'fold')
    
    % First Convolution
    %convolution2dLayer([10 1], 8, 'Name', 'conv1', 'NumChannels', 1, 'Stride', 10, 'Bias', zeros(1, 1, 8), 'BiasLearnRateFactor', 0)
    %batchNormalizationLayer('Name', 'bn1', 'Offset', zeros(1, 1, 8), 'OffsetLearnRateFactor', 0)
    %reluLayer('Name', 'relu1')
    %maxPooling2dLayer([5 1], 'Name', 'mpling1', 'Stride', 5)

    % Second Convolution
    %convolution2dLayer([5 1], 64, 'Name', 'conv2', 'NumChannels', 32, 'Stride', 3)
    %batchNormalizationLayer('Name', 'bn2')
    %reluLayer('Name', 'relu2')
    %maxPooling2dLayer([2 1], 'Name', 'mpling2', 'Stride', 1)
    
    % Third Convolution
    %convolution2dLayer([3 1], 128, 'Name', 'conv3', 'NumChannels', 64, 'Stride', 1)
    %batchNormalizationLayer('Name', 'bn3')
    %tanhLayer('Name', 'relu3')
    %maxPooling2dLayer([2 1], 'Name', 'mpling3', 'Stride', 1)
    
    % Flatten Layer
    %flattenLayer('Name', 'flatten1')
    
    % Fully Connected Layer
    fullyConnectedLayer(numHiddenUnits1, 'Name', 'fc2', 'Bias', zeros(50, 1), 'BiasLearnRateFactor', 0)
    dropoutLayer(0.1, 'Name', 'dpo1')
    reluLayer('Name', 'relu4')
    
    % First Upsampling
    %transposedConv2dLayer([3 1], 128, 'Name', 'tconv1', 'Stride', 1)
    
    % Second Upsampling
    %transposedConv2dLayer([10 1], 8, 'Name', 'tconv2', 'Stride', 10, 'Bias', zeros(1, 1, 8), 'BiasLearnRateFactor', 0)
    %batchNormalizationLayer('Name', 'bn4', 'Offset', zeros(1, 1, 8), 'OffsetLearnRateFactor', 0)
    %reluLayer('Name', 'relu3')
    
    % Third Upsampling
    %maxUnpooling2dLayer('Name', 'munpling')
    %transposedConv2dLayer([5 1], 32, 'Name', 'tconv3', 'Stride', 5)
    %reluLayer('Name', 'relu4')
    
    sequenceUnfoldingLayer('Name', 'unfold')
    
    % Second Flatten Layer
    flattenLayer('Name', 'flatten2')
    
    fullyConnectedLayer(numResponses3, 'Name', 'fc3', 'Bias', zeros(numResponses3, 1), 'BiasLearnRateFactor', 0)
    %dropoutLayer(0.1, 'Name', 'dpo1')
    
    regressionLayer('Name', 'rl')];



maxEpochs = 100;
miniBatchSize = 256;

%sLayer1 = skipLayer(128,4,64,'1');
%sLayer2 = skipLayer(256,4,128,'2');

lgraph = layerGraph(layers);
lgraph = connectLayers(lgraph,'fold/miniBatchSize','unfold/miniBatchSize');

%{
lgraph = addLayers(lgraph,sLayer1);
lgraph = addLayers(lgraph,sLayer2);
lgraph = connectLayers(lgraph,'relu','resizingConv1');
lgraph = connectLayers(lgraph,'lrelu2','resizingConv2');

lgraph = connectLayers(lgraph,'skippingBn1','add1/in2');
lgraph = connectLayers(lgraph,'skippingBn2','add2/in2');
%}

options = trainingOptions('adam', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'InitialLearnRate',0.001, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.5, ...
    'LearnRateDropPeriod',10, ...
    'GradientThreshold',Inf, ...
    'L2Regularization',0.001, ...
    'Shuffle','every-epoch', ...
    'Plots','training-progress', ...
    'Verbose',1, ...
    'ExecutionEnvironment', 'gpu');

%'L2Regularization',0.1, ...

net4 = trainNetwork(cData1, cResponse1, lgraph, options);

function layers = convolutionalUnit(numF,stride,numC,tag)
layers = [ ...
    convolution2dLayer(5,numF,'NumChannels',numC,'Padding','same','Stride',stride,'Name',[tag,'conv1'])
    batchNormalizationLayer('Name',[tag,'BN1'])
    tanhLayer('Name',[tag,'elu1'])
    
    convolution2dLayer(5,numF,'NumChannels',numC,'Padding','same','Stride',stride,'Name',[tag,'conv2'])
    batchNormalizationLayer('Name',[tag,'BN2'])];
    
end

function skippingLayer = skipLayer(numF,stride,numC,tag)
skippingLayer = [ ...
    
    convolution2dLayer([1 1], numF, 'Name', ['resizingConv' tag], 'NumChannels', numC, 'Stride', stride)
    batchNormalizationLayer('Name', ['skippingBn', tag])];

end