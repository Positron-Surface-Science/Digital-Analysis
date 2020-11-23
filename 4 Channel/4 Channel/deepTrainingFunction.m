numFeatures = 1951;
numHiddenUnits1 = 1000;
numHiddenUnits2 = 100;
numResponses1 = 500;
numResponses2 = 250;
numResponses3 = 1951;

layers = [ ...
    %imageInputLayer([438 1 1], 'Name', 'im')
    sequenceInputLayer([numFeatures 1 1], 'Name', 'input')
    %
    sequenceFoldingLayer('Name', 'fold')
    
    
    
    convolution2dLayer([5 1], 64, 'Name', 'conv', 'NumChannels', 1, 'Stride', 5)
    batchNormalizationLayer('Name', 'bn')
    tanhLayer('Name', 'relu')

    

    convolutionalUnit(128, 2, 64,'convolUnit')
    
    additionLayer(2, 'Name', 'add1')
    tanhLayer('Name', 'lrelu2')
    
    convolutionalUnit(256, 2, 128,'convolUnit2')

    additionLayer(2, 'Name', 'add2')
    tanhLayer('Name', 'lrelu3')


    averagePooling2dLayer([3 1], 'Name', 'mpling', 'Stride', 1)
    
%{
    fullyConnectedLayer(numHiddenUnits2, 'Name', 'fc2')
    dropoutLayer(0.5, 'Name', 'dpo1')
    tanhLayer('Name', 'lrelu4')
  %}  
    sequenceUnfoldingLayer('Name', 'unfold')
    
    flattenLayer('Name', 'flatten')

    lstmLayer(numHiddenUnits1, 'Name', 'lstm1')

    fullyConnectedLayer(numResponses3, 'Name', 'fc3')
    dropoutLayer(0.5, 'Name', 'dpo2')
    %tanhLayer('Name', 'relu3')
    regressionLayer('Name', 'rl')];



maxEpochs = 100;
miniBatchSize = 256;

sLayer1 = skipLayer(128,4,64,'1');
sLayer2 = skipLayer(256,4,128,'2');

lgraph = layerGraph(layers);
lgraph = addLayers(lgraph,sLayer1);
lgraph = addLayers(lgraph,sLayer2);
lgraph = connectLayers(lgraph,'relu','resizingConv1');
lgraph = connectLayers(lgraph,'lrelu2','resizingConv2');
lgraph = connectLayers(lgraph,'fold/miniBatchSize','unfold/miniBatchSize');
lgraph = connectLayers(lgraph,'skippingBn1','add1/in2');
lgraph = connectLayers(lgraph,'skippingBn2','add2/in2');

options = trainingOptions('adam', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'InitialLearnRate',0.01, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.5, ...
    'LearnRateDropPeriod',5, ...
    'GradientThreshold',Inf, ...
    'L2Regularization',0.00001, ...
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