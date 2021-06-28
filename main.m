dt = 0.01;

neuron = OLM_SingleCompartment_Modified(4000,dt);

neuron.I_stim(:) = 4/100;
neuron.eulerStabilize(6000,0.1);

neuron.I_stim(1000/dt+1:1000/dt+2000/dt) = neuron.I_stim(1000/dt+1:1000/dt+2000/dt)-120/100;
prg = waitbar(0,'Running Simulation');
while neuron.eulerStep() % Time stepper returns false once t == nsteps
    if (mod(neuron.t,round(neuron.nsteps/25))==0)
        waitbar(neuron.t/neuron.nsteps,prg);
    end
end
close(prg);
plot(neuron.V);
% hold on
% plot(1:length(neuron.V),neuron.V);
% NEURON = fopen('0.03pA.txt','r');
% V = fscanf(NEURON,'%f');
% plot(V);
% fclose(NEURON);
% hold off

