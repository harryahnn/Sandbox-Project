[sampleStruct, probStruct, Comments] = scfread('sample.scf');
figure
hold on
plot(sampleStruct.A);
plot(sampleStruct.C);
plot(sampleStruct.G);
plot(sampleStruct.T);
legend('A','C','G','T');