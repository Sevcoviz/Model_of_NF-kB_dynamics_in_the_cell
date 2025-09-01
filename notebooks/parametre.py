#parametre.py

class Parameters:
    def __init__(self):
        self.d_K = 0.15     # rýchlostná konštanta degradácie IKKα
        self.d = 3.0        # rýchlostná konštanta disociácie (rozpadu) komplexu NF-κB:IκBα
        self.gamma = 0.2    # rýchlostná konštanta pre spontánnu degradáciu komplexu NF-κB:IκBα
        self.d_I = 0.24     # rýchlostná konštanta spontánnej degradácie voľného IκBα
        self.p = 9.0        # faktor aktivity IKKα v rámci IKKα-indukovanej degradácie IκBα
        self.A = 200.0      # rýchlostná konštanta asociácie (viazania) voľného NF-κB a voľného IκBα
        self.kappa = 0.2    # škálovací faktor pre IKKα-indukovanú degradáciu voľného IκBα
        self.K_P = 16.0     # rýchlostná konštanta translácie IκBα mRNA na IκBα proteín
        self.d_R = 2.7      # rýchlostná konštanta degradácie IκBα mRNA
        self.k_ON = 7.5     # rýchlostná konštanta aktivácie génu IκBα transkripcie jadrovým NF-κB
        self.k_OFF = 15.0   # rýchlostná konštanta inaktivácie génu IκBα transkripcie jadrovým IκBα

        self.k_imp_N = 0.1  # rýchlostná konštanta importu voľného cytoplazmatického NF-κB do jadra
        self.n_assoc = 500  # rýchlostná konštanta asociácie (viazania) komplexu NF-κB:IκBα v jadre


# A single class to hold all parameters for clarity and easy modification.
# class CompartmentalParameters:
#     def __init__(self):
#         # Basal system parameters
#         self.d_K = 0.15      # IKKα degradation
#         self.d = 3.0         # {NI} complex dissociation
#         self.gamma = 0.2     # Spontaneous {NI} degradation
#         self.d_I = 0.24      # Free IκBα spontaneous degradation
#         self.p = 9.0         # IKKα-induced IκBα degradation factor
#         self.k_A = 200.0     # Association of NF-κB and IκBα in cytoplasm
#         self.k_A_N = 200.0   # Association of NF-κB and IκBα in nucleus
#         self.kappa = 0.2     # IKKα-induced degradation scaling factor
#         self.K_P = 16.0      # IκBα mRNA translation rate
#         self.d_R = 2.7       # IκBα mRNA degradation
#         self.k_ON = 7.5      # IκBα gene activation by nuclear NF-κB
#         self.k_OFF = 15.0    # IκBα gene inactivation by nuclear IκBα

#         # New compartmentalization parameters
#         self.k_imp_N = 0.5   # Free NF-κB import to nucleus
#         self.k_exp_NI = 0.01  # {NI} complex export from nucleus
#         self.k_imp_I = 0.28  # Free IκBα import to nucleus
#         self.k_exp_I = 0.28  # Free IκBα export from nucleus


class CompartmentalParameters:
    def __init__(self):
        # IKK dynamics
        self.k_IKK = 0.05       # IKKα synthesis (mM/s) - new parameter to sustain oscillations
        self.d_K = 0.005        # IKKα degradation (s^-1)

        # NF-κB:IκBα complex dynamics
        self.gamma = 0.0001     # spontaneous NI degradation (s^-1)
        self.d = 0.0001         # NFκB:IκBα complex dissociation (s^-1)

        # Spontaneous degradation rates
        self.d_I = 0.0001      # Free IκBα spontaneous degradation (s^-1)
        
        # IKK-induced degradation / catalytic strength
        self.p = 0.015          # IKKα-induced NI catalytic factor (mM^-1 s^-1)
        self.p_I = 0.005        # IKKα-induced free IκBα catalytic factor (mM^-1 s^-1) - accounts for direct degradation

        # Association (binding) rates
        self.k_A = 0.001          # cytoplasmic NFκB-IκB association (mM^-1 s^-1)
        self.k_A_N = 0.001        # nuclear NFκB-IκB association (mM^-1 s^-1)

        # mRNA / translation / degradation
        self.K_P = 0.0005       # IκBα mRNA -> protein translation rate (s^-1)
        self.d_R = 0.000025     # IκBα mRNA degradation (s^-1)

        # Gene switching rates (nucleus)
        self.k_ON = 0.0015      # gene activation by N_N (mM^-1 s^-1)
        self.k_OFF = 0.00025    # gene inactivation by I_N (mM^-1 s^-1)

        # Compartment transport
        self.k_imp_N = 0.00015    # NFκB nuclear import (s^-1)
        self.k_exp_N = 0.000001   # NFκB nuclear export (s^-1) - very slow for free NF-kB
        self.k_imp_I = 0.000001   # IκBα nuclear import (s^-1) - very slow for free IκB
        self.k_exp_I = 0.0001     # IκBα nuclear export (s^-1)
        self.k_imp_NI = 0.000001  # NI complex nuclear import (s^-1) - very slow for complex
        self.k_exp_NI = 0.0001    # NI complex nuclear export (s^-1)
