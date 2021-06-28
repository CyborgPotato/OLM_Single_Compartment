classdef OLM_SingleCompartment_Modified < handle
%OLM_SINGLECOMPARTMENT_MODIFIED Oriens-Lacunosum/Moleculare Single
%Compartment model
    properties
        %% Time dependent variables
        V % Membrane voltage
        I_NaT
        m_NaT
        h_NaT
        I_Kdrf
        m_Kdrf
        h_Kdrf 
        I_KA
        m_KA
        h_KA
        I_M
        m_M
        I_H
        r_H
        I_L
        I_stim % modified after creation for arbitrary stimulus
        t % Current timestep
        dt
        nsteps % number of timesteps
        %%Set Conductance values, units: pS/μm^2
        G_NaT = 70.986184915201491;
        G_Kdrf = 115.46932938891074;
        G_KA = 76.077776610698493;
        G_M = 0.13738940328219354;
        G_H = 0.106309;
        G_L = 0.075833;
        G_pas = 0.075833;
        %% Set Reversal potentials etc.
        E_Na = 90; % mV
        E_K = -95; % mV
        E_L = -64.6; % mV
        E_pas = -64.6; % mV
        E_h = -34.0056; % mV
        C_m = 0.27; % μF/cm^2
        T = 34; % Celsius
        qt = 3^((34-23)/10); % 34 == T
        SA = 31403.36016528357; % μm^2
        %% I_NA Constants
        V_Shift = -4.830346371483079; %mV
        %% I_Kdrf Constants
        tau_h = 1000;
        %% I_H Constants
        t1 = 8.5657797;
        t2 = 0.0296317;
        t3 = -6.9145;
        t4 = 0.1803;
        t5 = 4.3566601e-5;
        k = 9.9995804; %milliseconds
        V_12 = -103.69; %mV
    end
    
    methods (Access=protected)
        %% I_NaT related ODEs
        function I = Na(p,V,m,h)
            I = p.G_NaT*(m^3)*h*(V-p.E_Na);
        end
        function m = mdt_NaT(p,V,m1)
            am = (-0.1*(V+38-p.V_Shift)) / ...
                 (exp(-(V+38-p.V_Shift)/10)-1);
            bm = 4*exp(-(V+63-p.V_Shift)/18);
            
            m = am*(1-m1)-bm*m1;
        end
        function h = hdt_NaT(p,V,h1)
            ah = 0.07*exp(-(V+63-p.V_Shift)/20);
            bh = 1/(1+exp(-(V+33-p.V_Shift)/10));
            
            h = ah*(1-h1)-bh*h1;
        end
        %% I_Kdrf related ODEs
        function I = Kdrf(p,V,m,h)
            I = p.G_Kdrf*m*h*(V-p.E_K);
        end
        function m = mdt_Kdrf(p,V,m1)
            minf = (1/(1+exp(-(V+36.2)/16.1)))^4;
            taum = (27.8*exp((V+33)/14.3)) / ...
                   (p.qt*(1+exp((V+33)/10)));
            
            m = (minf-m1)/taum;
        end
        function h = hdt_Kdrf(p,V,h1)
            hinf = ((0.92)/(1+exp((V+40.6)/7.8))) + 0.08;
            tauh = p.tau_h; % 1000
            
            h = (hinf-h1)/tauh;
        end
        %% I_KA related ODEs
        function I = KA(p,V,m,h)
            I = p.G_KA*m*h*(V-p.E_K);
        end
        function m = mdt_KA(p,V,m1)
            minf = (1/(1+exp(-(V+41.4)/26.6)))^4;
            taum = 0.5/p.qt;
            
            m = (minf - m1)/taum;
        end
        function h = hdt_KA(p,V,h1)
            hinf = 1/(1+exp((V+78.5)/6));
            tauh = (0.17*(V+105))/(p.qt);
            
            h = (hinf-h1)/tauh;
        end
        %% I_M Related ODEs
        function I = M(p,V,m)
            I = p.G_M*m*(V-p.E_K);
        end
        function m = mdt_M(p,V,m1)
            minf = 1/(1+exp(-(V+27)/7));
            taum = 1 / ...
                      (0.003/exp(-(V+63)/15) + ...
                       0.003/exp((V+63)/15));
            
            m = (minf-m1)/taum;
        end
        %% I_H Related ODEs
        function I = H(p,V,r)
            I = p.G_H*r*(V-p.E_h);
        end
        function r = rdt_H(p,V,r1)
            rinf = 1/(1+exp((V-p.V_12)/(p.k)));
            tauh = (1/(exp(-p.t1-p.t2*V)+exp(-p.t3+p.t4*V))) + +p.t5;
            r = (rinf-r1)/(tauh);
        end
        %% I_L Current
        function I = L(p,V)
            I = p.G_L*(V-p.E_L);
        end
        %% Membrane Voltage ODE
        function v = Vdt(p,I_NaT,I_Kdrf,I_KA,I_M,I_H,I_L)
            v = (-I_NaT-I_Kdrf-I_KA-I_M-I_H-I_L+p.I_stim(p.t))/p.C_m;
        end
    end    
    
    methods (Static)
        function step = lin(x,dx,dt)
        % Takes a linear 'step' adjusting step = x+dx*dt;
            step = x+dx*dt;
        end
    end
    
    methods
        function p = OLM_SingleCompartment_Modified(tstop,dt)
        %OLM_SINGLECOMPARTMENT_MODIFIED constructor, based on modified
        %OLM model -- tstop and dt units are milliseconds
            %% Allocate tiem dependent variables
            p.nsteps = tstop/dt+1; % Number of time steps to calculate
            p.dt = dt;
            p.V = zeros(p.nsteps,1);
            p.I_NaT = zeros(p.nsteps,1);
            p.m_NaT = zeros(p.nsteps,1);
            p.h_NaT = zeros(p.nsteps,1);
            p.I_Kdrf = zeros(p.nsteps,1);
            p.m_Kdrf = zeros(p.nsteps,1);
            p.h_Kdrf = zeros(p.nsteps,1);
            p.I_KA = zeros(p.nsteps,1);
            p.m_KA = zeros(p.nsteps,1);
            p.h_KA = zeros(p.nsteps,1);
            p.I_M = zeros(p.nsteps,1);
            p.m_M = zeros(p.nsteps,1);
            p.I_H = zeros(p.nsteps,1);
            p.r_H = zeros(p.nsteps,1);
            p.I_L = zeros(p.nsteps,1);
            p.I_stim = zeros(p.nsteps,1);
            p.t=1;
        end
        
        function running = eulerStep(p)
        % Numerically intergrates one time step via Forward Euler
        % Use function in a while loop, allowing for interaction with
        % the model as it evolves, such as applying a Kalman filter
        % dynamically to a model
            
            running = true; % whether we have reached tstop/nsteps
            warning('off','');
            t = p.t; %#ok<*PROP> % ease of reference to current timestep
            dt = p.dt; % ease of reference to stepsize
            
            V = p.V(p.t); % ease of reference to current voltage
            
            % Step through each current before stepping membrane voltage
            %% I_NaT ODE
            p.m_NaT(t+1) = p.lin(p.m_NaT(t),mdt_NaT(p,V,p.m_NaT(t)),dt);
            p.h_NaT(t+1) = p.lin(p.h_NaT(t),hdt_NaT(p,V,p.h_NaT(t)),dt);
            p.I_NaT(t) = Na(p,V,p.m_NaT(t),p.h_NaT(t));
            %% I_Kdrf ODE
            p.m_Kdrf(t+1) = p.lin(p.m_Kdrf(t),mdt_Kdrf(p,V,p.m_Kdrf(t)),dt);
            p.h_Kdrf(t+1) = p.lin(p.h_Kdrf(t),hdt_Kdrf(p,V,p.h_Kdrf(t)),dt);
            p.I_Kdrf(t) = Kdrf(p,V,p.m_Kdrf(t),p.h_Kdrf(t));
            %% I_KA ODE
            p.m_KA(t+1) = p.lin(p.m_KA(t),mdt_KA(p,V,p.m_KA(t)),dt);
            p.h_KA(t+1) = p.lin(p.h_KA(t),hdt_KA(p,V,p.h_KA(t)),dt);
            p.I_KA(t) = KA(p,V,p.m_KA(t),p.h_KA(t));
            %% I_M ODE
            p.m_M(t+1) = p.lin(p.m_M(t),mdt_M(p,V,p.m_M(t)),dt);
            p.I_M(t) = M(p,V,p.m_M(t));
            %% I_H ODE
            p.r_H(t+1) = p.lin(p.r_H(t),rdt_H(p,V,p.r_H(t)),dt);
            p.I_H(t) = H(p,V,p.r_H(t));
            %% I_L
            p.I_L(t+1) = L(p,V);
            %% Membrane Voltage ODE
            p.V(t+1) = p.lin(V,Vdt(p,p.I_NaT(t),p.I_Kdrf(t),p.I_KA(t),p.I_M(t),p.I_H(t),p.I_L(t)),dt);
            p.t = p.t+1;
            if (p.t == p.nsteps)
                running = false;
            end
        end
        
        function running = rk4Step(p)
            running = true;
            
            %% Temporary variables to make calculation better
            V1 = p.V(p.t);
            V2 = p.V(p.t+1);
            %% Calculate k1 terms of all ODEs
            % I_NaT ODE
             % dm/dt
            k1_mNaT = mdt_NaT(p,p.V(p.t),p.m_NaT(p.t))*p.dt;
             % dh/dt
            k1_hNaT = hdt_NaT(p,p.V(p.t),p.h_NaT(p.t))*p.dt;
            % I_Kdrf ODE
             % dm/dt
            k1_mKdrf = mdt_Kdrf(p,p.V(p.t),p.m_Kdrf(p.t))*p.dt;
             % dh/dt
            k1_hKdrf = hdt_Kdrf(p,p.V(p.t),p.h_Kdrf(p.t))*p.dt;
            % I_KA ODE
             % dm/dt
            k1_mKA = mdt_KA(p,p.V(p.t),p.m_KA(p.t))*p.dt;
             % dh/dt
            k1_hKA = hdt_KA(p,p.V(p.t),p.h_KA(p.t))*p.dt;
            % I_M ODE
             % dm/dt
            k1_mM = mdt_M(p,p.V(p.t),p.m_M(p.t))*p.dt;
            % I_H ODE
             % dr/dt
            k1_rH = rdt_H(p,p.V(p.t),p.r_H(p.t))*p.dt;
            % Voltage Membrane ODE
            k1_V = Vdt(p)*p.dt/2;
            V2 = V1 + k1_V/3;
            V1 = V1 + k1_V;
            %% Calculate k2 terms of all ODEs
            % I_NaT ODE
             % dm/dt
            k2_mNaT = mdt_NaT(p,p.V(p.t)+0.5*k1_V,p.m_NaT(p.t)+k1_mNaT);
             % dh/dt
            k2_hNaT = hdt_NaT(p,p.V(p.t)+0.5*k1_V,p.h_NaT(p.t)+k1_hNaT);
            % I_Kdrf ODE
             % dm/dt
            k2_mKdrf = mdt_Kdrf(p,p.V(p.t)+0.5*k1_V,p.m_Kdrf(p.t)+k1_mKdrf);
             % dh/dt
            k2_hKdrf = hdt_Kdrf(p,p.V(p.t)+0.5*k1_V,p.h_Kdrf(p.t)+k1_hKdrf);
            % I_KA ODE
             % dm/dt
            k2_mKA = mdt_KA(p,p.V(p.t)+0.5*k1_V,p.m_KA(p.t)+k1_mKA);
             % dh/dt
            k2_hKA = hdt_KA(p,p.V(p.t)+0.5*k1_V,p.h_KA(p.t)+k1_hKA);
            % I_M ODE
             % dm/dt
            k2_mM = mdt_M(p,p.V(p.t)+0.5*k1_V,p.m_M(p.t)+k1_mM);
            % I_H ODE
             % dr/dt
            k2_rH = rdt_H(p,p.V(p.t)+0.5*k1_V,p.r_H(p.t)+k1_rH);
            % Voltage Membrane ODE
            k2_V = Vdt(p);
            
            %% Forward time step
            p.t = p.t+1;
            if (p.t == p.nsteps)
                running = false;
            end
        end
        
        function [] = eulerStabilize(p,tstop,dt)
            s = OLM_SingleCompartment_Modified(tstop,dt);
            s.V(1)=-74;
            s.I_stim(:) = p.I_stim(1);
            %% Copy all properties of p to s
            %%Set Conductance values, units: pS/μm^2
            s.G_NaT = p.G_NaT;
            s.G_Kdrf = p.G_Kdrf;
            s.G_KA = p.G_KA;
            s.G_M = p.G_M;
            s.G_H = p.G_H;
            s.G_L = p.G_L;
            s.G_pas = p.G_pas;
            %% Set Reversal potentials etc.
            s.E_Na = p.E_Na;
            s.E_K = p.E_K;
            s.E_L = p.E_L;
            s.E_pas = p.E_pas;
            s.E_h = p.E_h;
            s.C_m = p.C_m;
            s.T = p.T;
            s.qt = p.qt;
            s.SA = p.SA;
            %% I_NA Constants
            s.V_Shift = p.V_Shift;
            %% I_Kdrf Constants
            s.tau_h = p.tau_h;
            %% I_H Constants
            s.t1 = p.t1;
            s.t2 = p.t2;
            s.t3 = p.t3;
            s.t4 = p.t4;
            s.t5 = p.t5;
            s.k = p.k;
            s.V_12 = p.V_12;
            %% Stabilize currents
            while eulerStep(s)
            end
            %% Copy stablized currents and voltages to p
            p.V(1) = s.V(end);
            p.I_NaT(1) = s.I_NaT(end);
            p.m_NaT(1) = s.m_NaT(end);
            p.h_NaT(1) = s.h_NaT(end);
            p.I_Kdrf(1) = s.I_Kdrf(end);
            p.m_Kdrf(1) = s.m_Kdrf(end);
            p.h_Kdrf(1) = s.h_Kdrf(end);
            p.I_KA(1) = s.I_KA(end);
            p.m_KA(1) = s.m_KA(end);
            p.h_KA(1) = s.h_KA(end);
            p.I_M(1) = s.I_M(end);
            p.m_M(1) = s.m_M(end);
            p.I_H(1) = s.I_H(end);
            p.r_H(1) = s.r_H(end);
            p.I_L(1) = s.I_L(end);
            clear s
        end
    end
end

