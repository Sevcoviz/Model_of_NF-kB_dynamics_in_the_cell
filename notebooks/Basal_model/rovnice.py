from parametre import Parameters


class BasalSystem:
    def __init__(self, params=None):
        self.parameters = params if params else Parameters()

    def dK(self, K):
        return -self.parameters.d6 * K
    
    def dN(self, N, K, I):
        return self.parameters.a2 * (1 - N) + \
               self.parameters.d4 * self.parameters.d3 * (1 - N) + \
               self.parameters.d1 * K * (1 - N) - self.parameters.a1 * N * I

    def dI(self, N, K, I, R):
        return self.parameters.a2 * (1 - N) - \
               self.parameters.d2 * self.parameters.d1 * K * I - \
               self.parameters.a1 * N * I + \
               self.parameters.t4 * R - self.parameters.d3 * I

    def dR(self, R, G): #mRNA IkB
        return self.parameters.t3 * G - self.parameters.d5 * R

    def dG(self, N, G, I):
        return self.parameters.t1 * N * (1 - G) - self.parameters.t2 * I * G