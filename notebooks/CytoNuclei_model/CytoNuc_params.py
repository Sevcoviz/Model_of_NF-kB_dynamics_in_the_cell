class CytoNucParamsExact:
    def __init__(self, IKK_stimulation=0.5):

        self.a1 = 30.0       
        self.a2 = 0.03       
        self.a3 = 30.0       
        self.a4 = 0.03       

        self.IKK = IKK_stimulation
        self.d1 = 1.05 * self.IKK
        self.d5 = 0.017

        self.t3 = 1.03
        self.t4 = 0.24

        self.k1 = 5.4        
        self.k2 = 0.018      
        self.k3 = 0.012      
        self.k4 = 0.83      


