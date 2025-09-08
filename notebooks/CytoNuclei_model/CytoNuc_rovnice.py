
class NFkBSystemExact:
    """NF-κB signaling dynamics (exact article equations)."""

    def __init__(self, params):
        self.p = params

    def rhs(self, t, y):
        N, Nn, I, In, Im, NI, NIn = y
        p = self.p

        dNn = p.k1 * N - p.a3 * Nn * In + p.a4 * NIn
        dIm = p.t3 * (Nn ** 2) - p.d5 * Im
        dI = p.t4 * Im - p.a1 * N * I + p.a2 * NI - p.k2 * I + p.k3 * In
        dN = -p.a1 * N * I + (p.a2 + p.d1) * NI - p.k1 * N
        dNI = p.a1 * N * I - (p.a2 + p.d1) * NI + p.k4 * NIn
        dNIn = p.a3 * Nn * In - (p.a4 + p.k4) * NIn
        dIn = p.k2 * I - p.k3 * In - p.a3 * Nn * In + p.a4 * NIn

        return [dN, dNn, dI, dIn, dIm, dNI, dNIn]

