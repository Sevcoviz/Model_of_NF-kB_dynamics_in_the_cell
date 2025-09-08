class CytoNucParamsExact:
    def __init__(self, IKK_stimulation=0.5):
        # UNIFIED PARAMETER NAMES

        # Association / Dissociation
        self.a1 = 30.0       # Original: k_f - Binding / dissociation (μM^-1 min^-1)
        self.a2 = 0.03       # Original: k_b - Binding / dissociation (min^-1)
        self.a3 = 30.0       # Original: k_fn - Binding / dissociation (μM^-1 min^-1)
        self.a4 = 0.03       # Original: k_bn - Binding / dissociation (min^-1)

        # Degradation
        # IKK is now the variable stimulus parameter
        self.IKK = IKK_stimulation
        self.d1 = 1.05 * self.IKK  # Original: alpha - IKK-dependent degradation (min^-1)
        self.d5 = 0.017      # Original: gamma_m - mRNA degradation (min^-1)
        
        # Transcription / translation
        self.t3 = 1.03       # Original: k_t - Transcription / translation (μM^-1 min^-1)
        self.t4 = 0.24       # Original: k_tl - Transcription / translation (min^-1)

        # Transport
        self.k1 = 5.4        # Original: k_in - Transport (min^-1)
        self.k2 = 0.018      # Original: k_Iin - Transport (min^-1)
        self.k3 = 0.012      # Original: k_Iout - Transport (min^-1)
        self.k4 = 0.83       # Original: k_NIout - Transport (min^-1)


