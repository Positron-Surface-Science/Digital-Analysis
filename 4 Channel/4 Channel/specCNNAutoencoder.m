numFeatures1 = 1000;
numFeatures2 = 1;
numHiddenUnits1 = 500;
numHiddenUnits2 = 500;
numHiddenUnits3 = 1000;
numResponses = 1000;

layers = [ ...
    
    sequenceInputLayer([numFeatures1 numFeatures2 1], 'Name', 'input')
    %lstmLayer(numHiddenUnits2, 'Name', 'bll1')
    sequenceFoldingLayer('Name', 'fold')
    
    % First Convolution
    convolution2dLayer([4 1], 32, 'Name', 'conv1', 'NumChannels', 1, 'Stride', [1 1], 'Padding', [0 0])
    %batchNormalizationLayer('Name', 'bn1')
    leakyReluLayer('Name', 'relu1')
    %maxPooling2dLayer([2 1], 'Name', 'mpling1', 'Stride', 1)

    % Second Convolution
    convolution2dLayer([4 1], 64, 'Name', 'conv2', 'NumChannels', 32, 'Stride', [1 1], 'Padding', [0 0])
    %batchNormalizationLayer('Name', 'bn2')
    leakyReluLayer('Name', 'relu2')
    maxPooling2dLayer([2 1], 'Name', 'mpling1', 'Stride', [2 1], 'HasUnpoolingOutputs', true)
    
    % Third Convolution
    convolution2dLayer([4 1], 128, 'Name', 'conv3', 'NumChannels', 64, 'Stride', [1 1], 'Padding', [0 0])
    %batchNormalizationLayer('Name', 'bn2')
    leakyReluLayer('Name', 'relu3')
    maxPooling2dLayer([2 1], 'Name', 'mpling2', 'Stride', [2 1], 'HasUnpoolingOutputs', true)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CENTER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Second Upsampling
    %transposedConv2dLayer([50 1], 128, 'Name', 'tconv1', 'Stride', 1)
    maxUnpooling2dLayer('Name', 'munpling2')
    convolution2dLayer([4 1], 64, 'Name', 'conv4', 'NumChannels', 128, 'Stride', [1 1], 'Padding', [3 0])
    %dropoutLayer(0.5, 'Name', 'dpo1')
    %batchNormalizationLayer('Name', 'bn4')
    leakyReluLayer('Name', 'relu4')
    
    %additionLayer(2, 'Name', 'add2')
    
    % Third Upsampling
    %transposedConv2dLayer([122 1], 64, 'Name', 'tconv2', 'Stride', 1)
    maxUnpooling2dLayer('Name', 'munpling1')
    convolution2dLayer([4 1], 32, 'Name', 'conv5', 'NumChannels', 64, 'Stride', [1 1], 'Padding', [3 0])
    %dropoutLayer(0.5, 'Name', 'dpo2')
    leakyReluLayer('Name', 'relu5')
    %batchNormalizationLayer('Name', 'bn5')
    %)
    
    %additionLayer(2, 'Name', 'add1')
    
    %transposedConv2dLayer([249 1], 32, 'Name', 'tconv3', 'Stride', 1)
    convolution2dLayer([4 1], 1, 'Name', 'conv6', 'NumChannels', 32, 'Stride', [1 1], 'Padding', [0 0])
    %dropoutLayer(0.5, 'Name', 'dpo3')
    leakyReluLayer('Name', 'relu6')
    
    %fullyConnectedLayer(numHiddenUnits1, 'Name', 'fc1')
    %leakyReluLayer('Name', 'relu7')
    
    %fullyConnectedLayer(numHiddenUnits2, 'Name', 'fc2')
    %leakyReluLayer('Name', 'relu8')
    
    sequenceUnfoldingLayer('Name', 'unfold')
    %leakyReluLayer('Name', 'relu5')
    % Second Flatten Layer
    flattenLayer('Name', 'flatten2')
    

    fullyConnectedLayer(numResponses, 'Name', 'fc0')
    %leakyReluLayer('Name', 'tanh1')
    %dropoutLayer(0.1, 'Name', 'dpo2')
    
    regressionLayer('Name', 'rl')];

maxEpochs = 2000;
miniBatchSize = 128;

%sLayer1 = skipLayer(32,2,32,'1',[42 0]);
%sLayer2 = skipLayer(64,2,64,'2',[27 0]);

lgraph = layerGraph(layers);
lgraph = connectLayers(lgraph, 'fold/miniBatchSize', 'unfold/miniBatchSize');


%lgraph = addLayers(lgraph,sLayer1);
%lgraph = addLayers(lgraph,sLayer2);
%lgraph = connectLayers(lgraph,'relu1','resizingConv1');
%lgraph = connectLayers(lgraph,'relu2','resizingConv2');
%lgraph = connectLayers(lgraph, 'relu4', 'relu5');
%lgraph = connectLayers(lgraph,'resizingConv1','add1/in2');
%lgraph = connectLayers(lgraph,'resizingConv2','add2/in2');
%}
lgraph = connectLayers(lgraph,'mpling1/indices','munpling1/indices');
lgraph = connectLayers(lgraph,'mpling1/size','munpling1/size');

lgraph = connectLayers(lgraph,'mpling2/indices','munpling2/indices');
lgraph = connectLayers(lgraph,'mpling2/size','munpling2/size');

options = trainingOptions('adam', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'InitialLearnRate',0.0005, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.75, ...
    'LearnRateDropPeriod',1, ...
    'GradientThreshold',Inf, ...
    'L2Regularization',0.00001, ...
    'Shuffle','every-epoch', ...
    'Plots','training-progress', ...
    'Verbose',1, ...
    'ExecutionEnvironment', 'gpu');

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

function skippingLayer = skipLayer(numF, stride, numC, tag, pad)
skippingLayer = [ ...
    
    convolution2dLayer([1 1], numF, 'Name', ['resizingConv' tag], 'NumChannels', numC, 'Stride', stride, 'Padding', pad)];
    %batchNormalizationLayer('Name', ['skippingBn', tag])];

end
