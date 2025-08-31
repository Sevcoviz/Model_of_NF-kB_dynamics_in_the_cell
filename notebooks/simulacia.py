#simulacia.py
import numpy as np
from rovnice import BasalSystem
from parametre import Parameters

class Simulate_Basal_System:
    def __init__(self, K0=1.0, N0=0.0, I0=0.65, R0=0.0, G0=0.0, dt=0.01, t_max=5.0):
        self.K0 = K0
        self.N0 = N0
        self.I0 = I0
        self.R0 = R0
        self.G0 = G0
        self.dt = dt
        self.t_max = t_max

        self.system = BasalSystem()
        self.params = self.system.parameters
        self.n_steps = int(t_max / dt)

        # Outputs initialized to None
        self.t = None
        self.K = None
        self.N = None
        self.I = None
        self.R = None
        self.G = None
        self.G_mm = None

    def run(self):
        # Allocate arrays
        self.K = np.zeros(self.n_steps)
        self.N = np.zeros(self.n_steps)
        self.I = np.zeros(self.n_steps)
        self.R = np.zeros(self.n_steps)
        self.G = np.zeros(self.n_steps)
        self.G_mm = np.zeros(self.n_steps)

        # Initial values
        self.K[0] = self.K0
        self.N[0] = self.N0
        self.I[0] = self.I0
        self.R[0] = self.R0
        self.G[0] = self.G0

        self.G_mm[0] = (self.params.k_ON * self.N0) / (self.params.k_ON * self.N0 + self.params.k_OFF * self.I0)

        # Simulation loop
        for i in range(1, self.n_steps):
            dK_val = self.system.dK(self.K[i-1])
            dN_val = self.system.dN(self.N[i-1], self.K[i-1], self.I[i-1])
            dI_val = self.system.dI(self.N[i-1], self.K[i-1], self.I[i-1], self.R[i-1])
            dR_val = self.system.dR(self.R[i-1], self.G[i-1])
            dG_val = self.system.dG(self.N[i-1], self.G[i-1], self.I[i-1])

            self.G_mm[i] = (self.params.k_ON * self.N[i-1]) / (self.params.k_ON * self.N[i-1] + self.params.k_OFF * self.I[i-1])

            self.K[i] = self.K[i-1] + dK_val * self.dt
            self.N[i] = self.N[i-1] + dN_val * self.dt
            self.I[i] = self.I[i-1] + dI_val * self.dt
            self.R[i] = self.R[i-1] + dR_val * self.dt
            self.G[i] = self.G[i-1] + dG_val * self.dt

        self.t = np.linspace(0, self.t_max, self.n_steps)

    def get_results(self):
        return self.t, self.K, self.N, self.I, self.R, self.G, self.G_mm
    

