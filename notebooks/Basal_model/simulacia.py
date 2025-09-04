# #simulacia.py
# import numpy as np
# from rovnice import BasalSystem
# from parametre import Parameters
# from scipy.integrate import solve_ivp

# class Simulate_Basal_System:
#     def __init__(self, K0=1.0, N0=0.0, I0=0.65, R0=0.0, G0=0.0, dt=0.01, t_max=5.0):
#         self.K0 = K0
#         self.N0 = N0
#         self.I0 = I0
#         self.R0 = R0
#         self.G0 = G0
#         self.dt = dt
#         self.t_max = t_max

#         self.system = BasalSystem()
#         self.params = self.system.parameters
#         self.n_steps = int(t_max / dt)

#         # Outputs initialized to None
#         self.t = None
#         self.K = None
#         self.N = None
#         self.I = None
#         self.R = None
#         self.G = None
#         self.G_mm = None

#     def run(self):
#         # Allocate arrays
#         self.K = np.zeros(self.n_steps)
#         self.N = np.zeros(self.n_steps)
#         self.I = np.zeros(self.n_steps)
#         self.R = np.zeros(self.n_steps)
#         self.G = np.zeros(self.n_steps)
#         self.G_mm = np.zeros(self.n_steps)

#         # Initial values
#         self.K[0] = self.K0
#         self.N[0] = self.N0
#         self.I[0] = self.I0
#         self.R[0] = self.R0
#         self.G[0] = self.G0

#         self.G_mm[0] = (self.params.k_ON * self.N0) / (self.params.k_ON * self.N0 + self.params.k_OFF * self.I0)

#         # Simulation loop
#         for i in range(1, self.n_steps):
#             dK_val = self.system.dK(self.K[i-1])
#             dN_val = self.system.dN(self.N[i-1], self.K[i-1], self.I[i-1])
#             dI_val = self.system.dI(self.N[i-1], self.K[i-1], self.I[i-1], self.R[i-1])
#             dR_val = self.system.dR(self.R[i-1], self.G[i-1])
#             dG_val = self.system.dG(self.N[i-1], self.G[i-1], self.I[i-1])

#             self.G_mm[i] = (self.params.k_ON * self.N[i-1]) / (self.params.k_ON * self.N[i-1] + self.params.k_OFF * self.I[i-1])

#             self.K[i] = self.K[i-1] + dK_val * self.dt
#             self.N[i] = self.N[i-1] + dN_val * self.dt
#             self.I[i] = self.I[i-1] + dI_val * self.dt
#             self.R[i] = self.R[i-1] + dR_val * self.dt
#             self.G[i] = self.G[i-1] + dG_val * self.dt

#         self.t = np.linspace(0, self.t_max, self.n_steps)

#     def get_results(self):
#         return self.t, self.K, self.N, self.I, self.R, self.G, self.G_mm
    


import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from rovnice import BasalSystem
from parametre import Parameters


    # This class takes a system instance and provides methods to simulate
    # and plot the dynamics using the scipy ODE solver.
class BasalSystemSimulator:
    def __init__(self, system: BasalSystem):
        """
        Initializes the simulator with a BasalSystem instance.
        
        Args:
            system (BasalSystem): An instance of the BasalSystem class containing the model equations.
        """
        self.system = system
        self.params = system.parameters

    def _rhs(self, t, y):
        """
        Defines the right-hand side of the system of ordinary differential equations.
        This is the function that solve_ivp will call.
        
        Args:
            t (float): Time (not used in this autonomous system but required by solve_ivp).
            y (list or np.ndarray): The current state vector [K, N, I, R, G].
        
        Returns:
            list: The derivatives [dK/dt, dN/dt, dI/dt, dR/dt, dG/dt].
        """
        K, N, I, R, G = y

        dK_dt = self.system.dK(K)
        dN_dt = self.system.dN(N, K, I)
        dI_dt = self.system.dI(N, K, I, R)
        dR_dt = self.system.dR(R, G)
        dG_dt = self.system.dG(N, G, I)
        
        return [dK_dt, dN_dt, dI_dt, dR_dt, dG_dt]

    def simulate(self, t_span=(0, 5.0), y0=None, t_eval=None):
        """
        Runs the simulation using scipy's solve_ivp.
        
        Args:
            t_span (tuple, optional): The (start, end) time for the simulation. Defaults to (0, 5.0).
            y0 (list or np.ndarray, optional): Initial conditions for [K, N, I, R, G]. 
                                            If None, defaults are used.
            t_eval (np.ndarray, optional): Time points where the solution is stored. 
                                        If None, solve_ivp determines them.
        
        Returns:
            OdeResult: The solution object from solve_ivp.
        """
        if y0 is None:
            # Default initial conditions from the original notebook
            y0 = [1.0, 0.0, 0.65, 0.0, 0.0]  # K0, N0, I0, R0, G0
        
        if t_eval is None:
            t_eval = np.linspace(t_span[0], t_span[1], 500)

        sol = solve_ivp(
            fun=self._rhs,
            t_span=t_span,
            y0=y0,
            t_eval=t_eval,
            method="LSODA"  # A good general-purpose solver
        )
        return sol

    def plot_results(self, sol):
        """
        Plots the simulation results in a 2x2 grid.
        The first plot shows all variable dynamics together.
        
        Args:
            sol (OdeResult): The solution object returned by the simulate method.
        """
        # Extract variables from the solution object
        t = sol.t
        K, N, I, R, G = sol.y
        
        # Calculate the quasi-steady-state approximation for G
        G_mm = (self.params.t1 * N) / (self.params.t1 * N + self.params.t2 * I)

        # Create the plot grid
        fig, axs = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('Basal NF-κB System Dynamics', fontsize=16)

        # --- MODIFIED: Plot 1 now shows all dynamics ---
        axs[0, 0].plot(t, K, label='K (Active IKK)')
        axs[0, 0].plot(t, N, label='N (Nuclear NF-κB)')
        axs[0, 0].plot(t, I, label='I (Free IκBα)')
        axs[0, 0].plot(t, R, label='R (IκBα mRNA)')
        axs[0, 0].plot(t, G, label='G (Active Gene)')
        axs[0, 0].plot(t, G_mm, label='G_mm (QSS Approx.)', linestyle='--')
        axs[0, 0].set_title('Full System Dynamics')
        axs[0, 0].set_xlabel('Time')
        axs[0, 0].set_ylabel('Concentration / Fraction')
        axs[0, 0].legend(fontsize=8)
        axs[0, 0].grid(True)

        # Plot 2: Focused view on IκBα and its mRNA
        axs[0, 1].plot(t, I, label='I (Free IκBα)')
        axs[0, 1].plot(t, R, label='R (IκBα mRNA)')
        axs[0, 1].set_title('IκBα and its mRNA')
        axs[0, 1].set_xlabel('Time')
        axs[0, 1].set_ylabel('Concentration')
        axs[0, 1].legend()
        axs[0, 1].grid(True)

        # Plot 3: Focused view on Active IKK (K) and Nuclear NF-κB (N)
        axs[1, 0].plot(t, K, label='K (Active IKK)')
        axs[1, 0].plot(t, N, label='N (Nuclear NF-κB)')
        axs[1, 0].set_title('Active IKK and Nuclear NF-κB')
        axs[1, 0].set_xlabel('Time')
        axs[1, 0].set_ylabel('Concentration')
        axs[1, 0].legend()
        axs[1, 0].grid(True)

        # Plot 4: Focused view on Active IκBα Gene (G and G_mm)
        axs[1, 1].plot(t, G, label='G (Simulated Active Gene)')
        axs[1, 1].plot(t, G_mm, label='G_mm (Quasi-Steady-State)', linestyle='--')
        axs[1, 1].set_title('Active IκBα Gene Fraction')
        axs[1, 1].set_xlabel('Time')
        axs[1, 1].set_ylabel('Fraction')
        axs[1, 1].legend()
        axs[1, 1].grid(True)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.show()