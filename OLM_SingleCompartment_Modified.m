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
        %% Conudctances
        G_NaT
        G_Kdrf
        G_KA
        G_M
        G_H
        G_L
        G_pas
        %% Reversal Potentials etc.
        E_Na
        E_K
        E_L
        E_pas
        E_h
        C_m
        T
        qt
        SA % Surface Area
        %% I_NA Constants
        V_Shift
        %% I_Kdrf Constants
        tau_h
        %% I_H Constants
        t1
        t2
        t3
        t4
        t5
        k
        V_12
    end
    
    methods (Access=protected)
        %% I_NaT related ODEs
        function m = mdt_NaT(p,V,m1)
            am = -0.1*(V+38-p.V_Shift);
            am = am/(exp(-(V+38-p.V_Shift)/10)-1);
            bm = 4*exp(-(V+63-p.V_Shift)/18);
            
            m = am*(1-m1)-bm*m1;
        end
        function h = hdt_NaT(p,V,h1)
            ah = 0.07*exp(-(V+63-p.V_Shift)/20);
            bh = 1/( 1 + exp(-(V+33-p.V_Shift)/10));
            
            h = ah*(1-h1)-bh*h1;
        end
        %% I_Kdrf related ODEs
        function m = mdt_Kdrf(p,V,m1)
            minf = (1/(1+exp(-(V+36.2)/16.1)))^4;
            taum = 27.8*exp((V+33)/14.3);
            taum = taum/(p.qt*(1+exp((V+33)/10)));
            
            m = (minf-m1)/taum;
        end
        function h = hdt_Kdrf(p,V,h1)
            hinf = (0.92)/(1+exp((V+40.6)/7.8))+0.08;
            tauh = 1000;
            
            h = (hinf-h1)/tauh;
        end
        %% I_KA related ODEs
        function m = mdt_KA(p,V,m1)
            minf = (1/(1+exp(-(V+41.4)/26.6)))^4;
            taum = 0.5/p.qt;
            
            m = (minf-m1)/taum;
        end
        function h = hdt_KA(p,V,h1)
            hinf = 1/(1+exp((V+78.5)/6));
            tauh = 0.17*(V+105)/p.qt;
            
            h = (hinf-h1)/tauh;
        end
        %% I_M Related ODEs
        function m = mdt_M(p,V,m1)
            minf = 1/(1+exp(-(V+27)/7));
            taum = 1/(0.003*(1/exp(-(V+63)/15)+1/...
                             exp((V+63)/15)));
            
            m = (minf-m1)/taum;
        end
        %% I_H Related ODEs
        function r = rdt_H(p,V,r1)
            rinf = 1/(1+exp((V-p.V_12)/p.k));
            tauh = 1/(exp(-p.t1-p.t2*V)+...
                      exp(-p.t3+p.t4*V))+p.t5;
            
            r = (rinf-r1)/tauh;
        end
        %% Membrane Voltage ODE
        function v = Vdt(p)
            v = (-p.I_NaT(p.t)-p.I_Kdrf(p.t) - ...
                 p.I_KA(p.t)-p.I_M(p.t)-p.I_H(p.t) - ...
                 p.I_L(p.t)+p.I_stim(p.t))/p.C_m;
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
            %%Set Conductance values, units: pS/μm^2
            p.G_NaT = 70.986184915201491;
            p.G_Kdrf = 115.46932938891074;
            p.G_KA = 76.077776610698493;
            p.G_M = 0.13738940328219354;
            p.G_H = 0.106309;
            p.G_L = 0.075833;
            p.G_pas = 0.075833;
            %% Set Reversal potentials etc.
            p.E_Na = 90; % mV
            p.E_K = -96; % mV
            p.E_L = -64.6; % mV
            p.E_pas = -64.6; % mV
            p.E_h = -34.0056; % mV
            p.C_m = 0.27; % μF/cm^2
            p.T = 34; % Celsius
            p.qt = 3^((p.T-23)/10);
            p.SA = 31403.36016528357; % μm^2
            %% I_NA Constants
            p.V_Shift = -4.830346371483079; %mV
            %% I_Kdrf Constants
            p.tau_h = 1000;
            %% I_H Constants
            p.t1 = 8.5657797;
            p.t2 = 0.0296317;
            p.t3 = -6.9145;
            p.t4 = 0.1803;
            p.t5 = 4.3566601e-5;
            p.k = 9.9995804; %milliseconds
            p.V_12 = -103.69; %mV
        end
        
        function running = eulerStep(p)
        % Numerically intergrates one time step via Forward Euler
        % Use function in a while loop, allowing for interaction with
        % the model as it evolves, such as applying a Kalman filter
        % dynamically to a model
            running = true;
            
            % Step through each current before stepping membrane voltage
            %% I_NaT ODE
             % dm/dt
            p.m_NaT(p.t+1) = p.m_NaT(p.t)+mdt_NaT(p,p.V(p.t),p.m_NaT(p.t))*p.dt;
             % dh/dt
             p.h_NaT(p.t+1) = p.h_NaT(p.t)+hdt_NaT(p,p.V(p.t),p.h_NaT(p.t))*p.dt;
             % I_NaT
            p.I_NaT(p.t) = p.G_NaT*p.m_NaT(p.t)^3*p.h_NaT(p.t)*(p.V(p.t)-p.E_Na);
            %% I_Kdrf ODE
             % dm/dt
            p.m_Kdrf(p.t+1) = p.m_Kdrf(p.t)+mdt_Kdrf(p,p.V(p.t),p.m_Kdrf(p.t))*p.dt;
             % dh/dt
             p.h_Kdrf(p.t+1) = p.h_Kdrf(p.t)+hdt_Kdrf(p,p.V(p.t),p.h_Kdrf(p.t))*p.dt;
             % I_Kdrf
            p.I_Kdrf(p.t) = p.G_Kdrf*p.m_Kdrf(p.t)*p.h_Kdrf(p.t)*(p.V(p.t)-p.E_K);
            %% I_KA ODE
             % dm/dt
            p.m_KA(p.t+1) = p.m_KA(p.t)+ mdt_KA(p,p.V(p.t),p.m_KA(p.t))*p.dt;
             % dh/dt
             p.h_KA(p.t+1) = p.h_KA(p.t) + hdt_KA(p,p.V(p.t),p.h_KA(p.t))*p.dt;
             % I_KA
            p.I_KA(p.t) = p.G_KA*p.m_KA(p.t)*p.h_KA(p.t)*(p.V(p.t)-p.E_K);
            %% I_M ODE
             % dm/dt
            p.m_M(p.t+1) = p.m_M(p.t) + mdt_M(p,p.V(p.t),p.m_M(p.t))*p.dt;
             % I_M
            p.I_M(p.t) = p.G_M*p.m_M(p.t)*(p.V(p.t)-p.E_K);
            %% I_H ODE
             % dr/dt
            p.r_H(p.t+1) = p.r_H(p.t) + rdt_H(p,p.V(p.t),p.r_H(p.t))*p.dt;
             % I_H
            p.I_H(p.t) = p.G_H*p.r_H(p.t)*(p.V(p.t)-p.E_h);
            %% I_L
            p.I_L(p.t) = p.G_L*(p.V(p.t)-p.E_L);
            %% Membrane Voltage ODE
            p.V(p.t+1) = p.V(p.t) + Vdt(p)*p.dt;
            
            p.t = p.t+1;
            if (p.t == p.nsteps)
                running = false;
            end
        end
        
        function running = rk4Step(p)
            running = true;
            %% Calculate k1 terms of all ODEs
            % I_NaT ODE
             % dm/dt
            k1_mNaT = mdt_NaT(p,p.V(p.t),p.m_NaT(p.t));
             % dh/dt
            k1_hNaT = hdt_NaT(p,p.V(p.t),p.h_NaT(p.t));
            % I_Kdrf ODE
             % dm/dt
            k1_mKdrf = mdt_Kdrf(p,p.V(p.t),p.m_Kdrf(p.t));
             % dh/dt
            k1_hKdrf = hdt_Kdrf(p,p.V(p.t),p.h_Kdrf(p.t));
            % I_KA ODE
             % dm/dt
            k1_mKA = mdt_KA(p,p.V(p.t),p.m_KA(p.t));
             % dh/dt
            k1_hKA = hdt_KA(p,p.V(p.t),p.h_KA(p.t));
            % I_M ODE
             % dm/dt
            k1_mM = mdt_M(p,p.V(p.t),p.m_M(p.t));
            % I_H ODE
             % dr/dt
            k1_rH = rdt_H(p,p.V(p.t),p.r_H(p.t));
            % Voltage Membrane ODE
            %% TODO:Implement current dt to allow for step by step RK4
            k1_V = Vdt(p);
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
            s.V(1)=-80;
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
            figure
            plot(s.V)
            clear s
        end
    end
end

