#rovnice.py

from parametre import Parameters, CompartmentalParameters


class BasalSystem:
    def __init__(self, params=None):
        self.parameters = params if params else Parameters()

    def dK(self, K):
        return -self.parameters.d_K * K

    def dN(self, N, K, I):
        return self.parameters.d * (1 - N) + \
               self.parameters.gamma * self.parameters.d_I * (1 - N) + \
               self.parameters.p * K * (1 - N) - self.parameters.A * N * I

    def dI(self, N, K, I, R):
        return self.parameters.d * (1 - N) - \
               self.parameters.kappa * self.parameters.p * K * I - \
               self.parameters.A * N * I + \
               self.parameters.K_P * R - self.parameters.d_I * I

    def dR(self, R, G):
        return self.parameters.d_R * (G - R)

    def dG(self, N, G, I):
        return self.parameters.k_ON * N * (1 - G) - self.parameters.k_OFF * I * G
    


from parametre import CompartmentalParameters
# A single class to define the differential equations for the compartmentalized system.
# class CompartmentalSystem:
#     def __init__(self, params=None):
#         self.parameters = params if params else CompartmentalParameters()

#     def d_eqs(self, t, y):
#         """
#         Defines the system of ordinary differential equations (ODEs).
#         y is a vector of concentrations: [K, N_C, N_N, I_C, I_N, NI_C, NI_N, G, R]
#         """
#         K, N_C, N_N, I_C, I_N, NI_C, NI_N, G, R = y
#         p = self.parameters

#         # dK/dt
#         dK_dt = -p.d_K * K

#         # Rule: Free NF-κB is actively transported into the nucleus
#         # Rule: Free NF-κB is not transported out of the nucleus
#         dN_C_dt = p.k_exp_NI * NI_N + p.d * NI_C + p.p * K * NI_C - p.k_imp_N * N_C - p.k_A * N_C * I_C
#         dN_N_dt = p.k_imp_N * N_C - p.k_A_N * N_N * I_N

#         # Rule: IκB mRNA is translated to form IκB protein in the cytoplasm
#         # Rule: IκB protein can be transported in and out of the nucleus
#         dI_C_dt = p.K_P * R - p.k_A * N_C * I_C - p.k_imp_I * I_C + p.k_exp_I * I_N - p.d_I * I_C + p.p * K * NI_C + p.gamma * NI_C
#         dI_N_dt = -p.k_A_N * N_N * I_N + p.k_imp_I * I_C - p.k_exp_I * I_N - p.d_I * I_N

#         # Rule: In both compartments, IκB forms a complex with NF-κB
#         # Rule: {NI} complex cannot be imported into the nucleus
#         # Rule: {NI} complex can be exported out of the nucleus
#         dNI_C_dt = p.k_A * N_C * I_C - p.d * NI_C - p.p * K * NI_C - p.gamma * NI_C - p.k_exp_NI * NI_N
#         dNI_N_dt = p.k_A_N * N_N * I_N - p.k_exp_NI * NI_N

#         # Rule: NF-κB in the nucleus activates transcription of the IκBα gene
#         # The transcription rate is dependent on N_N.
#         dG_dt = p.k_ON * N_N * (1 - G) - p.k_OFF * I_N * G
#         dR_dt = p.d_R * (G - R)

#         return [dK_dt, dN_C_dt, dN_N_dt, dI_C_dt, dI_N_dt, dNI_C_dt, dNI_N_dt, dG_dt, dR_dt]



class CompartmentalSystem:
    def __init__(self, params=None):
        self.parameters = params if params else CompartmentalParameters()

    def ode_equations(self, t, y):
        """
        Defines the system of ODEs for the compartmentalized model.
        
        Args:
            t (float): Current time.
            y (list): A list of current concentrations [K, N_C, N_N, I_C, I_N, NI_C, NI_N, G, R].
        
        Returns:
            list: A list of the time derivatives of each concentration.
        """
        K, N_C, N_N, I_C, I_N, NI_C, NI_N, G, R = y
        p = self.parameters

        # dK/dt (IKKα dynamics)
        dK_dt = p.k_IKK - p.d_K * K

        # dN_C/dt (Cytoplasmic NF-κB)
        dN_C_dt = p.k_exp_NI * NI_N + p.d * NI_C + p.p * K * NI_C - p.k_imp_N * N_C - p.k_A * N_C * I_C

        # dN_N/dt (Nuclear NF-κB)
        dN_N_dt = p.k_imp_N * N_C - p.k_A_N * N_N * I_N

        # dI_C/dt (Cytoplasmic IκBα)
        dI_C_dt = p.K_P * R - p.k_A * N_C * I_C - p.k_imp_I * I_C + p.k_exp_I * I_N - p.d_I * I_C + p.p_I * K * I_C

        # dI_N/dt (Nuclear IκBα)
        dI_N_dt = -p.k_A_N * N_N * I_N + p.k_imp_I * I_C - p.k_exp_I * I_N + p.d_I * NI_N + p.d * NI_N

        # dNI_C/dt (Cytoplasmic Complex)
        dNI_C_dt = p.k_A * N_C * I_C - p.d * NI_C - p.p * K * NI_C - p.gamma * NI_C - p.k_imp_NI * NI_C + p.k_exp_NI * NI_N

        # dNI_N/dt (Nuclear Complex)
        dNI_N_dt = p.k_A_N * N_N * I_N + p.k_imp_NI * NI_C - p.k_exp_NI * NI_N - p.d * NI_N - p.d_I * NI_N

        # dG/dt (Gene State)
        dG_dt = p.k_ON * N_N * (1 - G) - p.k_OFF * I_N * G

        # dR/dt (IκBα mRNA)
        dR_dt = G - p.d_R * R

        return [dK_dt, dN_C_dt, dN_N_dt, dI_C_dt, dI_N_dt, dNI_C_dt, dNI_N_dt, dG_dt, dR_dt]