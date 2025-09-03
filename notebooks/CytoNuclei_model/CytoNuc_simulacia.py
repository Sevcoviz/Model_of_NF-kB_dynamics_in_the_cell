import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from CytoNuc_rovnice import NFkBSystemExact

# class NFkBSimulatorExact:
#     def __init__(self, system: NFkBSystemExact):
#         self.system = system

#     def simulate(self, t_span=(0, 500), y0=None, t_eval=None):
#         if y0 is None:
#             # Initial conditions: N=1 μM, IKK=0.5 μM, others zero
#             y0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#         if t_eval is None:
#             t_eval = np.linspace(t_span[0], t_span[1], 2000)

#         sol = solve_ivp(
#             fun=self.system.rhs,
#             t_span=t_span,
#             y0=y0,
#             t_eval=t_eval,
#             method="LSODA"
#         )
#         return sol

#     def plot(self, sol):
#         plt.figure(figsize=(10, 6))
#         plt.plot(sol.t, sol.y[1], label="Nn (nuclear NF-κB)", linewidth=2)
#         plt.plot(sol.t, sol.y[0], label="N (cytoplasmic NF-κB)")
#         plt.plot(sol.t, sol.y[2], label="I (cytoplasmic IκB)")
#         plt.plot(sol.t, sol.y[3], label="In (nuclear IκB)")
#         plt.plot(sol.t, sol.y[4], label="Im (mRNA IκB)")
#         plt.plot(sol.t, sol.y[5], label="NI (cytoplasmic complex)")
#         plt.plot(sol.t, sol.y[6], label="NIn (nuclear complex)")

#         plt.xlabel("Time (min)")
#         plt.ylabel("Concentration (μM)")
#         plt.legend()
#         plt.title("Exact NF-κB Signaling Dynamics (Article Model)")
#         plt.tight_layout()
#         plt.show()

#     def plot_totals(self, sol):
#         # Extract variables
#         N, Nn, I, In, Im, NI, NIn = sol.y

#         # Compute totals
#         NFkB_total = N + Nn
#         IkB_total  = I + In
#         NI_total   = NI + NIn
#         Im_total   = Im

#         # Plot
#         plt.figure(figsize=(10, 6))
#         plt.plot(sol.t, NFkB_total, label="Total NF-κB", linewidth=2)
#         plt.plot(sol.t, IkB_total, label="Total IκB", linewidth=2)
#         plt.plot(sol.t, NI_total, label="Total NF-κB/IκB Complex", linewidth=2)
#         plt.plot(sol.t, Im_total, label="IκB mRNA", linewidth=2, linestyle="--")

#         plt.xlabel("Time (min)")
#         plt.ylabel("Concentration (μM)")
#         plt.legend()
#         plt.title("Absolute Totals (NF-κB, IκB, Complex, mRNA)")
#         plt.tight_layout()
#         plt.show()

#     def plot_N_I_dynamics(self, sol):
#         # Extract variables
#         N, Nn, I, In, Im, NI, NIn = sol.y

#         # Compute totals
#         NFkB_total = N + Nn
#         IkB_total  = I + In
#         NI_total   = NI + NIn
#         Im_total   = Im

#         # Plot
#         plt.figure(figsize=(10, 6))
#         plt.plot(sol.t, NFkB_total, label="Total NF-κB", linewidth=2)
#         plt.plot(sol.t, IkB_total, label="Total IκB", linewidth=2)
#         # plt.plot(sol.t, NI_total, label="Total NF-κB/IκB Complex", linewidth=2)
#         # plt.plot(sol.t, Im_total, label="IκB mRNA", linewidth=2, linestyle="--")

#         plt.xlabel("Time (min)")
#         plt.ylabel("Concentration (μM)")
#         plt.legend()
#         plt.title("Absolute Totals (NF-κB, IκB)")

#         plt.tight_layout()
#         plt.show()

#     def plot_N_dynamics(self, sol):
#         # Extract variables
#         N, Nn, I, In, Im, NI, NIn = sol.y

#         # Plot
#         plt.figure(figsize=(10, 6))
#         plt.plot(sol.t, N, label="N (cytoplasmic NF-κB)", linewidth=2)
#         plt.plot(sol.t, Nn, label="Nn (nuclear NF-κB)", linewidth=2)
#         plt.xlabel("Time (min)")
#         plt.ylabel("Concentration (μM)")
#         plt.legend()
#         plt.title("NF-κB Dynamics (Cytoplasmic and Nuclear)")
#         plt.tight_layout()
#         plt.show()


class NFkBSimulatorExact:
    def __init__(self, system: NFkBSystemExact):
        self.system = system

    def simulate(self, t_span=(0, 500), y0=None, t_eval=None):
        if y0 is None:
            # Initial conditions: N=1 μM, IKK=0.5 μM, others zero
            y0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        if t_eval is None:
            t_eval = np.linspace(t_span[0], t_span[1], 2000)

        sol = solve_ivp(
            fun=self.system.rhs,
            t_span=t_span,
            y0=y0,
            t_eval=t_eval,
            method="LSODA"
        )
        return sol

    def plot_all(self, sol):
        # Extract variables
        N, Nn, I, In, Im, NI, NIn = sol.y

        # Totals
        NFkB_total = N + Nn
        IkB_total  = I + In
        NI_total   = NI + NIn
        Im_total   = Im

        # Create grid: 2 columns x 2 rows
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        axes = axes.flatten()

        # --- Panel 1: Detailed dynamics ---
        axes[0].plot(sol.t, Nn, label="Nn (nuclear NF-κB)", linewidth=2)
        axes[0].plot(sol.t, N, label="N (cytoplasmic NF-κB)")
        axes[0].plot(sol.t, I, label="I (cytoplasmic IκB)")
        axes[0].plot(sol.t, In, label="In (nuclear IκB)")
        axes[0].plot(sol.t, Im, label="Im (mRNA IκB)")
        axes[0].plot(sol.t, NI, label="NI (cytoplasmic complex)")
        axes[0].plot(sol.t, NIn, label="NIn (nuclear complex)")
        axes[0].set_title("Full Dynamics")
        axes[0].legend(fontsize=8)

        # --- Panel 2: Totals ---
        axes[1].plot(sol.t, NFkB_total, label="Total NF-κB", linewidth=2)
        axes[1].plot(sol.t, IkB_total, label="Total IκB", linewidth=2)
        axes[1].plot(sol.t, NI_total, label="Total Complex", linewidth=2)
        axes[1].plot(sol.t, Im_total, label="IκB mRNA", linestyle="--", linewidth=2)
        axes[1].set_title("Total Concentrations")
        axes[1].legend(fontsize=8)

        # --- Panel 3: NF-κB vs IκB totals ---
        axes[2].plot(sol.t, NFkB_total, label="Total NF-κB", linewidth=2)
        axes[2].plot(sol.t, IkB_total, label="Total IκB", linewidth=2)
        axes[2].set_title("NF-κB vs IκB (Totals)")
        axes[2].legend(fontsize=8)

        # --- Panel 4: NF-κB compartments ---
        axes[3].plot(sol.t, N, label="N (cytoplasmic)", linewidth=2)
        axes[3].plot(sol.t, Nn, label="Nn (nuclear)", linewidth=2)
        axes[3].set_title("NF-κB Cytoplasmic vs Nuclear")
        axes[3].legend(fontsize=8)

        # Format all axes
        for ax in axes:
            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Concentration (μM)")
            ax.grid(True)

        fig.suptitle("NF-κB Signaling Dynamics", fontsize=14, y=1.02)
        fig.tight_layout()
        plt.show()
