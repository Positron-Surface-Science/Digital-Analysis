numFeatures1 = 1000;
numFeatures2 = 1;
numHiddenUnits1 = 1000;
numHiddenUnits2 = 2000;
numHiddenUnits3 = 1000;
numResponses = 1000;
sizeFilters = [3 1];
sizePFilters = [2 1];
sizePStride = [2 1];

layers = [ ...
    
    sequenceInputLayer([numFeatures1 1 1], 'Name', 'input')
    sequenceFoldingLayer('Name', 'fold')
    
    % First Convolution
    convolution2dLayer([7 1], 128, 'Name', 'conv1', 'NumChannels', 1, 'Stride', [1 1], 'Padding', [3 0])
    %dropoutLayer(0.25, 'Name', 'dpo1')
    %batchNormalizationLayer('Name', 'bn1')
    leakyReluLayer(0.25, 'Name', 'relu1')
    
    convolution2dLayer(sizeFilters,128, 'Name', 'conv2', 'NumChannels', 128, 'Stride', [1 1], 'Padding', [1 0])
    %dropoutLayer(0.25, 'Name', 'dpo2')
    %batchNormalizationLayer('Name', 'bn2')
    leakyReluLayer(0.25, 'Name', 'relu2')
    maxPooling2dLayer(sizePFilters, 'Name', 'mpling1', 'Stride', sizePStride, 'HasUnpoolingOutputs', true)

    % Second Convolution
    convolution2dLayer(sizeFilters, 256, 'Name', 'conv3', 'NumChannels', 128, 'Stride', [1 1], 'Padding', [1 0])
    %dropoutLayer(0.25, 'Name', 'dpo3')
    %batchNormalizationLayer('Name', 'bn3')
    leakyReluLayer(0.25, 'Name', 'relu3')
    
    convolution2dLayer(sizeFilters, 256, 'Name', 'conv4', 'NumChannels', 256, 'Stride', [1 1], 'Padding', [1 0])
    %dropoutLayer(0.25, 'Name', 'dpo4')
    %batchNormalizationLayer('Name', 'bn4')
    leakyReluLayer(0.25, 'Name', 'relu4')
    maxPooling2dLayer(sizePFilters, 'Name', 'mpling2', 'Stride', sizePStride, 'HasUnpoolingOutputs', true)
    
    % Third Convolution
    convolution2dLayer(sizeFilters, 512, 'Name', 'conv5', 'NumChannels', 256, 'Stride', [1 1], 'Padding', [1 0])
    %dropoutLayer(0.25, 'Name', 'dpo5')
    %batchNormalizationLayer('Name', 'bn5')
    leakyReluLayer(0.25, 'Name', 'relu5')
    
    convolution2dLayer(sizeFilters, 512, 'Name', 'conv6', 'NumChannels', 512, 'Stride', [1 1], 'Padding', [1 0])
    %dropoutLayer(0.25, 'Name', 'dpo5')
    %batchNormalizationLayer('Name', 'bn5')
    leakyReluLayer(0.25, 'Name', 'relu6')
    maxPooling2dLayer(sizePFilters, 'Name', 'mpling3', 'Stride', sizePStride, 'HasUnpoolingOutputs', true)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CENTER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Second Upsampling
    maxUnpooling2dLayer('Name', 'munp3')
    
    convolution2dLayer(sizeFilters, 512, 'Name', 'conv7', 'NumChannels', 512, 'Stride', [1 1], 'Padding', [1 0])
    dropoutLayer(0.5, 'Name', 'dpo7')
    %batchNormalizationLayer('Name', 'bn6')
    leakyReluLayer(0.25, 'Name', 'relu7')
    
    convolution2dLayer(sizeFilters, 256, 'Name', 'conv8', 'NumChannels', 512, 'Stride', [1 1], 'Padding', [1 0])
    dropoutLayer(0.5, 'Name', 'dpo8')
    %batchNormalizationLayer('Name', 'bn6')
    leakyReluLayer(0.25, 'Name', 'relu8')
    
    % Third Upsampling
    maxUnpooling2dLayer('Name', 'munp2')
    
    additionLayer(2, 'Name', 'add2')
    
    convolution2dLayer(sizeFilters, 128, 'Name', 'conv9', 'NumChannels', 256, 'Stride', [1 1], 'Padding', [1 0])
    dropoutLayer(0.5, 'Name', 'dpo9')
    %batchNormalizationLayer('Name', 'bn7')
    leakyReluLayer(0.25,'Name', 'relu9')
    
    convolution2dLayer(sizeFilters, 128, 'Name', 'conv10', 'NumChannels', 128, 'Stride', [1 1], 'Padding', [1 0])
    dropoutLayer(0.5, 'Name', 'dpo10')
    %batchNormalizationLayer('Name', 'bn8')
    leakyReluLayer(0.1, 'Name', 'relu10')
    
    maxUnpooling2dLayer('Name', 'munp1')
    
    additionLayer(2, 'Name', 'add1')
    
    %convolution2dLayer([8 1], 32, 'Name', 'conv9', 'NumChannels', 32, 'Stride', [1 1], 'Padding', [3 0])
    %dropoutLayer(0.5, 'Name', 'dpo4')
    %batchNormalizationLayer('Name', 'bn6')
    %reluLayer('Name', 'relu9')
    
    convolution2dLayer(sizeFilters, 1, 'Name', 'conv11', 'NumChannels', 128, 'Stride', [1 1], 'Padding', [1 0])
    dropoutLayer(0.5, 'Name', 'dpo11')
    %batchNormalizationLayer('Name', 'bn9')
    leakyReluLayer(0.25, 'Name', 'relu11')
    
    crop2dLayer('centercrop', 'Name', 'crop1')
    
    additionLayer(2, 'Name', 'add0')
    
    fullyConnectedLayer(numHiddenUnits1, 'Name', 'fc1')
    dropoutLayer(0.5, 'Name', 'dpo12')
    leakyReluLayer(0.5, 'Name', 'relu12')
    
    fullyConnectedLayer(numHiddenUnits2, 'Name', 'fc2')
    %dropoutLayer(0.25, 'Name', 'dpo13')
    %reluLayer('Name', 'relu12')
    
    %globalAveragePooling2dLayer('Name','gapool1')
    
    sequenceUnfoldingLayer('Name', 'unfold')
    
    % Second Flatten Layer
    flattenLayer('Name', 'flatten2')

    fullyConnectedLayer(numResponses, 'Name', 'fc0')
    
    regressionLayer('Name', 'rl')];

maxEpochs = 2000;
miniBatchSize = 64;

%sLayer1 = skipLayer(64,1,64,'1',[0 0]);
%sLayer2 = skipLayer(128,1,128,'2',[0 0]);

lgraph = layerGraph(layers);
lgraph = connectLayers(lgraph, 'fold/miniBatchSize', 'unfold/miniBatchSize');


%lgraph = addLayers(lgraph,sLayer1);
%lgraph = addLayers(lgraph,sLayer2);
lgraph = connectLayers(lgraph,'conv2','add1/in2'); %'resizingConv1');
lgraph = connectLayers(lgraph,'conv4','add2/in2'); %resizingConv2');
%lgraph = connectLayers(lgraph, 'relu4', 'relu5');
%lgraph = connectLayers(lgraph,'resizingConv1','add1/in2');
%lgraph = connectLayers(lgraph,'resizingConv2','add2/in2');

lgraph = connectLayers(lgraph,'fold/out','add0/in2');

lgraph = connectLayers(lgraph,'fold/out','crop1/ref')
%}
lgraph = connectLayers(lgraph,'mpling1/indices','munp1/indices');
lgraph = connectLayers(lgraph,'mpling1/size','munp1/size');

lgraph = connectLayers(lgraph,'mpling2/indices','munp2/indices');
lgraph = connectLayers(lgraph,'mpling2/size','munp2/size');

lgraph = connectLayers(lgraph,'mpling3/indices','munp3/indices');
lgraph = connectLayers(lgraph,'mpling3/size','munp3/size');

options = trainingOptions('adam', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'InitialLearnRate',0.001, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.75, ...
    'LearnRateDropPeriod',2, ...
    'GradientThreshold',1, ...
    'L2Regularization',0.01, ...
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

function skippingLayer = skipLayer(numF, stride, numC, tag, pad)
skippingLayer = [ ...
    
    convolution2dLayer([1 1], numF, 'Name', ['resizingConv' tag], 'NumChannels', numC, 'Stride', stride, 'Padding', pad)];
    %batchNormalizationLayer('Name', ['skippingBn', tag])];

end
