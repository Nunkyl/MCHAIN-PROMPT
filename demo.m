function [] = demo()
%DEMO Model conformational movement for a Hha-H-NS46 charge zipper complex 
% We demonstrate routines of the MCHAIN-PROMPT package by modeling
% conformational motion between two conformations of a protein with 
% PDB ID 2MW2. To launch the example, use the following command: 
% *demo*.

file = 'samples/2MW2.pdb';
[initialPDBstruct, firstConf, lastConf] = getdata(file); 

nConf = 8;
model = trmcreate(firstConf, lastConf, nConf);

angleIndices = trmdistantangleindices(model);
f = @(x) trmobjfunc(model, angleIndices, x);

initial_point = createinitialpoint(model, angleIndices);

iterNum = 2;

options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'iter');
options = optimoptions(options,'MaxIter', iterNum);
options = optimoptions(options,'MaxFunEvals', Inf);
options = optimoptions(options,'GradObj','off');
options = optimoptions(options,'UseParallel',true);

[lowerBound, upperBound] = getbounds(model, angleIndices);

if max(size(gcp)) == 0 
    parpool
end

x = fmincon(f,initial_point,[],[],[],[],lowerBound, upperBound,[],...
    options);

modelOptimized = modelopt(model, angleIndices, x);

delete(gcp);

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

name = 'samples/2MW2_opt100.pdb';
toolboxPath = fileparts(which('demo'));
newPDB = trm2pdb(modelOptimized, initialPDBstruct);
pdbwrite(fullfile(toolboxPath, name), newPDB);

end
