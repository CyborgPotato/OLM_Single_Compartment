neuron = OLM_SingleCompartment_Modified(4000,0.1);

neuron.I_stim(:) = 4/neuron.SA;

neuron.eulerStabilize(6000,0.1);

neuron.I_stim(1000/0.1+1:1000/0.1+2000/0.1) = neuron.I_stim(1000/0.1+1:1000/0.1+2000/0.1)-120/neuron.SA;
% neuron.I_stim(1000/0.1+1:1000/0.1+2000/0.1) = neuron.I_stim(1000/0.1+1:1000/0.1+2000/0.1)-120/100;
% neuron.I_stim(1000/0.01+1:1000/0.01+2000/0.01) = neuron.I_stim(1000/0.01+1:1000/0.01+2000/0.01)+60/100;

while neuron.eulerStep() % Time stepper returns false once t == nsteps
    
end

plot(1:neuron.nsteps,neuron.V)