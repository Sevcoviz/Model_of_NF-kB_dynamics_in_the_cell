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



    #compartmental
        self.kdeg = 0.125 # degradace IKK
        self.k_bind_I = 0.2   # IKK*IkBa association  Fit
        self.k_bind = 0.5     # NF-κB:IκBα association

        self.k_imp_I = 0.28 # rýchlostná konštanta importu voľného cytoplazmatického IκBα do jadra
        self.k_exp_I = 0.07 # rýchlostná konštanta exportu voľného jadrového IκBα z jadra do cytoplazmy

        self.k_export_complex = 0.02
        self.k_dissociate = 0.05
        self.k_trans = 1.0 # rýchlostná konštanta translácie IκBα mRNA na IκBα proteín
        self.k_transc = 0.5 # mRNA transcription rate
        self.k_mRNA_export = 0.1 # rýchlostná konštanta exportu IκBα mRNA z jadra do cytoplazmy

        # self.k_ON = 0.1
        # self.k_OFF = 0.05
        self.Vc = 5
        self.Vn = 1
