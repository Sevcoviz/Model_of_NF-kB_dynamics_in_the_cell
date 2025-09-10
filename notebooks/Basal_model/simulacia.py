import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from rovnice import BasalSystem
from parametre import Parameters


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
            y0 = [1.0, 0.0, 0.65, 0.0, 0.0]  # K0, N0, I0, R0, G0
        
        if t_eval is None:
            t_eval = np.linspace(t_span[0], t_span[1], 500)

        sol = solve_ivp(
            fun=self._rhs,
            t_span=t_span,
            y0=y0,
            t_eval=t_eval,
            method="LSODA"  
        )
        return sol

    def plot_results(self, sol):
        """
        Plots the simulation results in a 2x2 grid.
        The first plot shows all variable dynamics together.
        
        Args:
            sol (OdeResult): The solution object returned by the simulate method.
        """
        t = sol.t
        K, N, I, R, G = sol.y
        
        G_mm = (self.params.t1 * N) / (self.params.t1 * N + self.params.t2 * I)

        fig, axs = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('Basal NF-κB System Dynamics', fontsize=16)

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

        axs[0, 1].plot(t, I, label='I (Free IκBα)')
        axs[0, 1].plot(t, R, label='R (IκBα mRNA)')
        axs[0, 1].set_title('IκBα and its mRNA')
        axs[0, 1].set_xlabel('Time')
        axs[0, 1].set_ylabel('Concentration')
        axs[0, 1].legend()
        axs[0, 1].grid(True)

        axs[1, 0].plot(t, K, label='K (Active IKK)')
        axs[1, 0].plot(t, N, label='N (Nuclear NF-κB)')
        axs[1, 0].set_title('Active IKK and Nuclear NF-κB')
        axs[1, 0].set_xlabel('Time')
        axs[1, 0].set_ylabel('Concentration')
        axs[1, 0].legend()
        axs[1, 0].grid(True)

        axs[1, 1].plot(t, G, label='G (Simulated Active Gene)')
        axs[1, 1].plot(t, G_mm, label='G_mm (Quasi-Steady-State)', linestyle='--')
        axs[1, 1].set_title('Active IκBα Gene Fraction')
        axs[1, 1].set_xlabel('Time')
        axs[1, 1].set_ylabel('Fraction')
        axs[1, 1].legend()
        axs[1, 1].grid(True)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.show()