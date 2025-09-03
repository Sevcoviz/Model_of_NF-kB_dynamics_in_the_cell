class CytoNucParamsExact:
    def __init__(self):
        # Transport
        self.k_in = 5.4       # min^-1
        self.k_Iin = 0.018    # min^-1
        self.k_Iout = 0.012   # min^-1
        self.k_NIout = 0.83   # min^-1

        # Transcription / translation
        self.k_t = 1.03       # μM^-1 min^-1
        self.k_tl = 0.24      # min^-1

        # Binding / dissociation
        self.k_f = 30.0       # μM^-1 min^-1
        self.k_fn = 30.0      # μM^-1 min^-1
        self.k_b = 0.03       # min^-1
        self.k_bn = 0.03      # min^-1

        # IKK-dependent degradation
        self.IKK = 0.5        # μM initial
        self.alpha = 1.05 * self.IKK  # min^-1

        # mRNA degradation
        self.gamma_m = 0.017  # min^-1
