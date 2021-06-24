neuron = OLM_SingleCompartment_Modified(5000,0.1);

neuron.eulerStabilize(2000,0.1)

neuron.I_stim(round(end/5):round(2*end/5)) = 60/neuron.SA*100;

while neuron.eulerStep() % Time stepper returns false once t+1 == nsteps
    
end

plot(1:neuron.nsteps,neuron.V)