#rovnice.py

# from parametre import Parameters


# class BasalSystem:
#     def __init__(self, params=None):
#         self.parameters = params if params else Parameters()

#     def dK(self, K):
#         return -self.parameters.d_K * K
    
#     def dN(self, N, K, I):
#         return self.parameters.d * (1 - N) + \
#                self.parameters.gamma * self.parameters.d_I * (1 - N) + \
#                self.parameters.p * K * (1 - N) - self.parameters.A * N * I

#     def dI(self, N, K, I, R):
#         return self.parameters.d * (1 - N) - \
#                self.parameters.kappa * self.parameters.p * K * I - \
#                self.parameters.A * N * I + \
#                self.parameters.K_P * R - self.parameters.d_I * I

#     # def dR(self, R, G): #mRNA IkB
#     #     return self.parameters.d_R * (G - R)

#     def dR(self, R, G): #mRNA IkB
#         return self.parameters.k_R * G - self.parameters.d_R * R

#     def dG(self, N, G, I):
#         return self.parameters.k_ON * N * (1 - G) - self.parameters.k_OFF * I * G
    


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