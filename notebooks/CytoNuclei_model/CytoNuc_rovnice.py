class NFkBSystemExact:
    """NF-κB signaling dynamics (exact article equations)."""

    def __init__(self, params):
        self.p = params

    def rhs(self, t, y):
        N, Nn, I, In, Im, NI, NIn = y
        p = self.p

        dNn = p.k_in * N - p.k_fn * Nn * In + p.k_bn * NIn
        dIm = p.k_t * (Nn ** 2) - p.gamma_m * Im
        dI = p.k_tl * Im - p.k_f * N * I + p.k_b * NI - p.k_Iin * I + p.k_Iout * In
        dN = -p.k_f * N * I + (p.k_b + p.alpha) * NI - p.k_in * N
        dNI = p.k_f * N * I - (p.k_b + p.alpha) * NI + p.k_NIout * NIn
        dNIn = p.k_fn * Nn * In - (p.k_bn + p.k_NIout) * NIn
        dIn = p.k_Iin * I - p.k_Iout * In - p.k_fn * Nn * In + p.k_bn * NIn

        return [dN, dNn, dI, dIn, dIm, dNI, dNIn]

