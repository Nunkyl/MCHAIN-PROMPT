function [] = demo2(file,name)

file0 = 'samples/insulin/new_insulin.pdb'; 
file1 = 'samples/calmodulin/calmodulin.pdb';
file2 = 'samples/short_insulin.pdb';
file3 = 'samples/2MP2/2MP2.pdb';
file4 = 'samples/2MXR.pdb';
file5 = 'samples/2LE9.pdb';
file6 = 'samples/KcsA.pdb';
file7 = 'samples/hemoglobin.pdb';
file8 = 'samples/final_KcsA.pdb';
file9 = 'samples/2LWW/2LWW.pdb';
file10 = 'samples/2RVB/2RVB.pdb';
file11 = 'samples/2myj/2myj.pdb';
file12 = 'samples/2m3x/2m3x.pdb';
file13 = 'samples/2m6i/2m6i.pdb';
file14 = 'samples/2k1e/2k1e.pdb';
file15 = 'samples/2n1f/2n1f.pdb';
file16 = 'samples/2mxu/2mxu.pdb';
file17 = 'samples/2mse/2mse.pdb';
file18 = 'samples/2mse/2mse_clean.pdb';
file19 = 'samples/2n1f/2n1f_altered.pdb';
file20 = 'samples/2mse/2mse_normal.pdb';


[initialPDBstruct, firstConf, lastConf] = getData(file);
nConf = 8;

separatechains = cell(length(firstConf),1);

for i=1:length(firstConf)
    separatechains{i} = optimize(trmcreate(firstConf{i}, lastConf{i}, ...
        nConf), 1000, 'on');
end

model = trmcreate(firstConf, lastConf, nConf);

for i=1:length(firstConf)
    model(i).psi = separatechains{i}.psi;
end

modelOptimized = optimize(model, 300, 'off');

toolboxPath = fileparts(which('calmodulin_demo'));
newPDB = trm2pdb(modelOptimized, initialPDBstruct);
pdbwrite(fullfile(toolboxPath, name), newPDB);

end

function [modelOptimized] = optimize(model, iterNum, grad)

angleIndices = trmdistantangleindices(model);
f = @(x) trmobjfunc(model, angleIndices, x);

initial_point = createInitialPoint(model, angleIndices);

options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'iter');%off
options = optimoptions(options,'MaxIter', iterNum);
options = optimoptions(options,'MaxFunEvals', Inf);
options = optimoptions(options,'GradObj',grad);
options = optimoptions(options,'UseParallel',true);

[lowerBound, upperBound] = getBounds(model, angleIndices);

if max(size(gcp)) == 0 
    parpool
end

x = fmincon(f,initial_point,[],[],[],[],lowerBound, upperBound,[],...
    options);

modelOptimized = modelOpt(model, angleIndices, x);

delete(gcp);

end

function [] = plotangles(model)
    trmplottranglediff(model);
    nAngles = ones(length(model))*50;
    trmplottranglediff(model, true, nAngles, '-o');
end

function [] = plotrmsds(model, modelOptimized)
 
    figure;
    trmplotadjrmsd({model, modelOptimized});
    legend('original', 'optimized');
    grid;
    title('RMSDs Between Adjacent Configurations');

    figure;
    confNo = 1;
    trmplotfixedrmsd({model, modelOptimized}, confNo);
    legend('original', 'optimized');
    grid;
    title('RMSDs to the First Configuration');

end

% angles = [140, 160];
% nAngles = getNumAng(model, angles);
% transformations = trmvariateinterpolations(model, nAngles);
