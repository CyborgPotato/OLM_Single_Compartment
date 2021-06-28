dt = 0.1;

neuron = OLM_SingleCompartment_Modified(4000,dt);

neuron.stim = 30/neuron.SA*1000;

RES = neuron.ode45Step();

% figure
hold on
plot(RES(:,1));
NEURON = fopen('0.03pA.txt','r');
V = fscanf(NEURON,'%f');
plot(V);
fclose(NEURON);
hold off

