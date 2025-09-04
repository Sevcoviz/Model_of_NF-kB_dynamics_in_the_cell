#parametre.py

# class Parameters:
#     def __init__(self):
#         self.d_K = 0.15     # rýchlostná konštanta degradácie IKKα
#         self.d = 3.0        # rýchlostná konštanta disociácie (rozpadu) komplexu NF-κB:IκBα
#         self.gamma = 0.2    # rýchlostná konštanta pre spontánnu degradáciu komplexu NF-κB:IκBα
#         self.d_I = 0.24     # rýchlostná konštanta spontánnej degradácie voľného IκBα
#         self.p = 9.0        # faktor aktivity IKKα v rámci IKKα-indukovanej degradácie IκBα
#         self.A = 200.0      # rýchlostná konštanta asociácie (viazania) voľného NF-κB a voľného IκBα
#         self.kappa = 0.2    # škálovací faktor pre IKKα-indukovanú degradáciu voľného IκBα
#         self.K_P = 16.0     # rýchlostná konštanta translácie IκBα mRNA na IκBα proteín
#         self.d_R = 2.7      # rýchlostná konštanta degradácie IκBα mRNA
#         self.k_ON = 7.5     # rýchlostná konštanta aktivácie génu IκBα transkripcie jadrovým NF-κB
#         self.k_OFF = 15.0   # rýchlostná konštanta inaktivácie génu IκBα transkripcie jadrovým IκBα
#         self.k_R = 2.5      # rýchlostná konštanta transkripcie génu IκBα na IκBα mRNA (nový parameter)

#         self.k_imp_N = 0.1  # rýchlostná konštanta importu voľného cytoplazmatického NF-κB do jadra
#         self.n_assoc = 500  # rýchlostná konštanta asociácie (viazania) komplexu NF-κB:IκBα v jadre




class Parameters:
    def __init__(self):
        # UNIFIED PARAMETER NAMES

        # Association / Dissociation
        self.a1 = 200.0      # Original: A - rýchlostná konštanta asociácie (viazania) voľného NF-κB a voľného IκBα
        self.a2 = 3.0        # Original: d - rýchlostná konštanta disociácie (rozpadu) komplexu NF-κB:IκBα

        # Degradation
        self.d1 = 9.0        # Original: p - faktor aktivity IKKα v rámci IKKα-indukovanej degradácie IκBα
        self.d2 = 0.2        # Original: kappa - škálovací faktor pre IKKα-indukovanú degradáciu voľného IκBα
        self.d3 = 0.24       # Original: d_I - rýchlostná konštanta spontánnej degradácie voľného IκBα
        self.d4 = 0.2        # Original: gamma - rýchlostná konštanta pre spontánnu degradáciu komplexu NF-κB:IκBα
        self.d5 = 2.7        # Original: d_R - rýchlostná konštanta degradácie IκBα mRNA
        self.d6 = 0.15       # Original: d_K - rýchlostná konštanta degradácie IKKα

        # Transcription / Translation
        self.t1 = 7.5        # Original: k_ON - rýchlostná konštanta aktivácie génu IκBα transkripcie jadrovým NF-κB
        self.t2 = 15.0       # Original: k_OFF - rýchlostná konštanta inaktivácie génu IκBα transkripcie jadrovým IκBα
        self.t3 = 2.5        # Original: k_R - rýchlostná konštanta transkripcie génu IκBα na IκBα mRNA
        self.t4 = 16.0       # Original: K_P - rýchlostná konštanta translácie IκBα mRNA na IκBα proteín