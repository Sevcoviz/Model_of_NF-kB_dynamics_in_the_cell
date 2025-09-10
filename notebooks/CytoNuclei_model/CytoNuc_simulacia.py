import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from CytoNuc_rovnice import NFkBSystemExact
from CytoNuc_params import CytoNucParamsExact

class NFkBSimulatorExact:
    def __init__(self, system: NFkBSystemExact):
        self.system = system

    def simulate(self, t_span=(0, 1000), y0=None, t_eval=None):
        if y0 is None:
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


    def plot_dynamics(self, sol):
        """Plots a comprehensive 4x2 dashboard of system dynamics."""
        N, Nn, I, In, Im, NI, NIn = sol.y
        
        # --- CALCULATIONS ---
        NFkB_total = N + Nn
        IkB_total = I + In
        NI_total = NI + NIn
        
        # Calculation for Fractional Activation plot
        NFkB_total_free = N + Nn
        fractional_Nn = Nn / (NFkB_total_free + 1e-9)

        # NEW: Calculations for Compartmental Pools and Complex plots
        NFkB_cyto_pool = N + NI
        NFkB_nuc_pool = Nn + NIn
        
        # --- PLOTTING SETUP ---
        # Changed to a 4x2 grid and increased vertical size
        fig, axes = plt.subplots(4, 2, figsize=(12, 16)) 
        axes = axes.flatten()

        # --- ROW 1: Original Plots ---
        # Panel 1 (axes[0]): Detailed dynamics
        axes[0].plot(sol.t, Nn, label="Nn (nuclear NF-κB)", linewidth=2)
        axes[0].plot(sol.t, N, label="N (cytoplasmic NF-κB)")
        axes[0].plot(sol.t, I, label="I (cytoplasmic IκB)")
        axes[0].plot(sol.t, In, label="In (nuclear IκB)")
        axes[0].plot(sol.t, Im, label="Im (mRNA IκB)")
        axes[0].plot(sol.t, NI, label="NI (cytoplasmic complex)")
        axes[0].plot(sol.t, NIn, label="NIn (nuclear complex)")
        axes[0].set_title("Full Dynamics")
        
        # Panel 2 (axes[1]): Totals
        axes[1].plot(sol.t, NFkB_total, label="Total NF-κB", linewidth=2)
        axes[1].plot(sol.t, IkB_total, label="Total IκB", linewidth=2)
        axes[1].plot(sol.t, NI_total, label="Total Complex", linewidth=2)
        axes[1].plot(sol.t, Im, label="IκB mRNA", linestyle="--", linewidth=2)
        axes[1].set_title("Total Concentrations")
        
        # --- ROW 2: Original Plots ---
        # Panel 3 (axes[2]): NF-κB vs IκB totals
        axes[2].plot(sol.t, NFkB_total, label="Total NF-κB", linewidth=2)
        axes[2].plot(sol.t, IkB_total, label="Total IκB", linewidth=2)
        axes[2].set_title("NF-κB vs IκB (Totals)")

        # Panel 4 (axes[3]): NF-κB compartments
        axes[3].plot(sol.t, N, label="N (cytoplasmic)", linewidth=2)
        axes[3].plot(sol.t, Nn, label="Nn (nuclear)", linewidth=2)
        axes[3].set_title("NF-κB Cytoplasmic vs Nuclear")

        # --- ROW 3: Advanced Plots from previous request ---
        # Panel 5 (axes[4]): Fractional Activation of NF-κB
        axes[4].plot(sol.t, fractional_Nn, color='orangered', linewidth=2, label='Fractional Nn')
        axes[4].set_title("Fractional Nuclear Activation of NF-κB")
        axes[4].set_ylabel("Fraction of Nn / (N + Nn)")
        axes[4].set_ylim(0, 1.05)

        # Panel 6 (axes[5]): Phase Portrait (Nn vs Im)
        axes[5].plot(Nn, Im, color='purple')
        axes[5].set_title("Phase Portrait: Activator vs. Feedback Source")
        axes[5].set_xlabel("Nuclear NF-κB (Nn) [μM]")
        axes[5].set_ylabel("IκB mRNA (Im) [μM]")

        # --- ROW 4: NEW Mechanistic Plots ---
        # Panel 7 (axes[6]): Total NF-κB Compartmental Pools
        axes[6].plot(sol.t, NFkB_cyto_pool, label="Total Cyto NF-κB (N+NI)", linewidth=2)
        axes[6].plot(sol.t, NFkB_nuc_pool, label="Total Nuc NF-κB (Nn+NIn)", linewidth=2)
        axes[6].set_title("Total NF-κB Compartmental Pools")
        axes[6].set_ylabel("Total Concentration (μM)")

        # Panel 8 (axes[7]): Dynamics of Complex Formation
        axes[7].plot(sol.t, NI, label="Cytoplasmic Complex (NI)", linewidth=2)
        axes[7].plot(sol.t, NIn, label="Nuclear Complex (NIn)", linewidth=2)
        axes[7].set_title("Dynamics of Complex Formation")
        axes[7].set_ylabel("Complex Concentration (μM)")
        
        # --- FORMATTING ---
        # Apply formatting carefully to avoid overwriting labels
        time_series_axes_indices = [0, 1, 2, 3, 4, 6, 7]
        concentration_axes_indices = [0, 1, 2, 3, 6, 7]
        
        for i in time_series_axes_indices:
            axes[i].set_xlabel("Time (min)")
            axes[i].grid(True)
            if i != 4: # Don't add legend to fractional plot if it only has one line
                axes[i].legend(fontsize=8)

        for i in concentration_axes_indices:
            axes[i].set_ylabel("Concentration (μM)")

        # Specific formatting for the non-time-series plot
        axes[5].grid(True)
                
        fig.suptitle(f"NF-κB Signaling Dynamics (IKK={self.system.p.IKK:.2f} µM)", fontsize=16, y=1.02)
        fig.tight_layout()
        plt.show()


    def _run_sensitivity_analysis_single_param(self, param_name, param_range, ikk_stim=0.5):
        """Generic internal method for running sensitivity analysis on one parameter."""
        peak_Nn_values = []
        auc_Nn_values = []
        final_Nn_values = []

        print(f"Running sensitivity analysis for {len(param_range)} '{param_name}' levels...")

        for val in param_range:
            params = CytoNucParamsExact(IKK_stimulation=ikk_stim)
            setattr(params, param_name, val)  
            self.system.p = params
            
            sol = self.simulate(t_span=(0, 1000))
            
            t, Nn = sol.t, sol.y[1]
            peak_Nn_values.append(np.max(Nn))
            auc_Nn_values.append(np.trapz(Nn, t))
            final_Nn_values.append(Nn[-1])
            
        print(f"Analysis for '{param_name}' finished.")
        return param_range, peak_Nn_values, auc_Nn_values, final_Nn_values

    def run_pathology_monitoring(self, IKK_stim_range):
        return self._run_sensitivity_analysis_single_param('IKK', IKK_stim_range, ikk_stim=None) # Special case

    def run_sensitivity_analysis_t3(self, t3_range, ikk_stim=0.5):
        return self._run_sensitivity_analysis_single_param('t3', t3_range, ikk_stim)
        
    def run_sensitivity_analysis_k1(self, k1_range, ikk_stim=0.5):
        """NEW: Runs sensitivity analysis for k1."""
        return self._run_sensitivity_analysis_single_param('k1', k1_range, ikk_stim)

    def run_sensitivity_analysis_k2(self, k2_range, ikk_stim=0.5):
        """NEW: Runs sensitivity analysis for k2."""
        return self._run_sensitivity_analysis_single_param('k2', k2_range, ikk_stim)




    def plot_advanced_dynamics(self, sol):
        """NEW: Plots more advanced visualizations of system dynamics."""
        N, Nn, I, In, Im, NI, NIn = sol.y
        
        NFkB_total_pool = N + Nn + NI + NIn  
        IkB_total_pool = I + In + NI + NIn     
        NFkB_total_free = N + Nn              
        fractional_Nn = Nn / (NFkB_total_free + 1e-9) 
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()

        axes[0].plot(Nn, Im, color='purple')
        axes[0].set_title("Phase Portrait: Activator vs. Feedback Source")
        axes[0].set_xlabel("Nuclear NF-κB (Nn) [μM]")
        axes[0].set_ylabel("IκB mRNA (Im) [μM]")
        axes[0].grid(True)
        
        axes[1].plot(Nn, In, color='green')
        axes[1].set_title("Phase Portrait: Activator vs. Nuclear Inhibitor")
        axes[1].set_xlabel("Nuclear NF-κB (Nn) [μM]")
        axes[1].set_ylabel("Nuclear IκB (In) [μM]")
        axes[1].grid(True)
        
        axes[2].plot(sol.t, fractional_Nn, color='orangered', linewidth=2)
        axes[2].set_title("Fractional Nuclear Activation of NF-κB")
        axes[2].set_xlabel("Time (min)")
        axes[2].set_ylabel("Fraction of Nn / (N + Nn)")
        axes[2].set_ylim(0, 1.05)
        axes[2].grid(True)

        axes[3].plot(sol.t, NFkB_total_pool, label="Total NF-κB Pool", linewidth=2)
        axes[3].plot(sol.t, IkB_total_pool, label="Total IκB Pool", linewidth=2)
        axes[3].set_title("Conservation of Total Protein")
        axes[3].set_xlabel("Time (min)")
        axes[3].set_ylabel("Total Concentration (μM)")
        axes[3].legend()
        axes[3].grid(True)
            
        fig.suptitle(f"Advanced NF-κB Dynamics (IKK={self.system.p.IKK:.2f} µM)", fontsize=16, y=1.02)
        fig.tight_layout()
        plt.show()
        
    def plot_pathology_report(self, IKK_stim_range, peak_Nn_values, auc_Nn_values, final_Nn_values):
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        axes[0].plot(IKK_stim_range, peak_Nn_values, 'o-', color='r')
        axes[0].set_title('Peak Nuclear NF-κB (Nn)')
        axes[0].set_ylabel('Peak Concentration (μM)')
        
        axes[1].plot(IKK_stim_range, auc_Nn_values, 'o-', color='g')
        axes[1].set_title('Sustained Nn Activation (AUC)')
        axes[1].set_ylabel('Total Activity (μM * min)')

        axes[2].plot(IKK_stim_range, final_Nn_values, 'o-', color='b')
        axes[2].set_title('Final Nn Steady-State')
        axes[2].set_ylabel('Final Concentration (μM)')

        for ax in axes:
            ax.grid(True)
            ax.set_xlabel('IKK Stimulus Concentration (μM)')
        
        fig.suptitle("Pathological Response Monitoring Pipeline", fontsize=16, y=1.02)
        fig.tight_layout()
        plt.show()

    def plot_sensitivity_report_t3(self, t3_range, peak_Nn_values, auc_Nn_values, final_Nn_values):
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        axes[0].plot(t3_range, peak_Nn_values, 'o-', color='r')
        axes[0].set_title('Peak Nuclear NF-κB (Nn)')
        axes[0].set_ylabel('Peak Concentration (μM)')
        
        axes[1].plot(t3_range, auc_Nn_values, 'o-', color='g')
        axes[1].set_title('Sustained Nn Activation (AUC)')
        axes[1].set_ylabel('Total Activity (μM * min)')

        axes[2].plot(t3_range, final_Nn_values, 'o-', color='b')
        axes[2].set_title('Final Nn Steady-State')
        axes[2].set_ylabel('Final Concentration (μM)')

        for ax in axes:
            ax.grid(True)
            ax.set_xlabel('IκB Transcription Rate (t3)')
        
        fig.suptitle("Sensitivity Analysis: Impact of Feedback Loop Strength (t3)", fontsize=16, y=1.02)
        fig.tight_layout()
        plt.show()
        
    def plot_sensitivity_report_transport(self, k1_results, k2_results):
        """NEW: Plots a comparative report for k1 and k2 sensitivity."""
        fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=False)
        
        k1_range, peak_k1, auc_k1, final_k1 = k1_results
        k2_range, peak_k2, auc_k2, final_k2 = k2_results
        
        k1_norm = k1_range / 5.4  # Baseline k1 is 5.4
        k2_norm = k2_range / 0.018 # Baseline k2 is 0.018
        
        axes[0].plot(k1_norm, peak_k1, 'o-', color='purple', label='NF-κB Import (k1)')
        axes[0].plot(k2_norm, peak_k2, 's-', color='orange', label='IκB Import (k2)')
        axes[0].set_title('Peak Nn Sensitivity')
        axes[0].set_ylabel('Peak Concentration (μM)')
        
        axes[1].plot(k1_norm, auc_k1, 'o-', color='purple', label='NF-κB Import (k1)')
        axes[1].plot(k2_norm, auc_k2, 's-', color='orange', label='IκB Import (k2)')
        axes[1].set_title('Sustained Activation (AUC) Sensitivity')
        axes[1].set_ylabel('Total Activity (μM * min)')

        axes[2].plot(k1_norm, final_k1, 'o-', color='purple', label='NF-κB Import (k1)')
        axes[2].plot(k2_norm, final_k2, 's-', color='orange', label='IκB Import (k2)')
        axes[2].set_title('Final Steady-State Sensitivity')
        axes[2].set_ylabel('Final Concentration (μM)')

        for ax in axes:
            ax.grid(True)
            ax.set_xlabel('Parameter Fold Change from Baseline')
            ax.legend()
        
        fig.suptitle("Sensitivity Analysis: Nuclear Import of Activator (k1) vs. Inhibitor (k2)", fontsize=16, y=1.02)
        fig.tight_layout()
        plt.show()



    def run_pathology_monitoring(self, IKK_stim_range):
        """Runs simulations across a range of IKK stimulus values and returns metrics."""
        peak_Nn_values = []
        auc_Nn_values = []
        final_Nn_values = []

        print(f"Running pathology pipeline for {len(IKK_stim_range)} stimulus levels...")

        for ikk_val in IKK_stim_range:
            self.system.p = CytoNucParamsExact(IKK_stimulation=ikk_val)
            sol = self.simulate(t_span=(0, 1000))
            
            t, Nn = sol.t, sol.y[1]
            peak_Nn_values.append(np.max(Nn))
            auc_Nn_values.append(np.trapz(Nn, t))
            final_Nn_values.append(Nn[-1])

        print("Pipeline finished.")
        return IKK_stim_range, peak_Nn_values, auc_Nn_values, final_Nn_values

        






