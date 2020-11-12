

import numpy as np

class Janaf():
  # Read Janaf Table from https://janaf.nist.gov/ as text
    def __init__(self, filename, MaxTemp = 1500):
        self.janaftable = open(filename, 'r')
        self.tabletexts = open(filename, 'r').readlines()
        self.tableText = open(filename, 'r').read()
        self.MaxTemp = MaxTemp
        self.PointTemp =np.array([])
        self.PointEntropy =np.array([])
        self.PointEnthalpy =np.array([])
        self.FitEnthalpy =np.array([])
        self.FitTempRange =np.array([])
        self.FitEntropy =np.array([])
        self.FitEntropyR2 =0
        self.FitEnthalpyR2 =0
        
    def get_datas(self):
        PointTemerature = [];  PointEntropy = []; PointEnthalpy = []; 
        for n,line in enumerate(self.tabletexts):
            temp = line.split()
            if n < 2: pass
            elif n == 2:
                ground_enthalpy= float(temp[4])  
            else:
                PointTemerature.append(float(temp[0]))
                PointEntropy.append(   float(temp[2]) * 0.0103642688  / 1000. )
                PointEnthalpy.append(( float(temp[4]) - ground_enthalpy )  * 0.0103642688 ) 
        self.PointTemp     = np.array(PointTemerature)
        self.PointEnthalpy = np.array(PointEnthalpy)
        self.PointEntropy  = np.array(PointEntropy)
        
    def fit_datas(self,FitPoints = 300):
        
        def r2(points, fitted):
            bar  = np.sum(points)/len(points)
            ssres = np.sum( (points- fitted)**2. )
            ssreg = np.sum( (fitted- bar)**2. )
            sstot = np.sum( (points - bar)**2. )
            return 1 - ssres/sstot
        
        self.FitTempRange = np.linspace(min(self.PointTemp), self.MaxTemp,FitPoints )
        
        self.fEntropyPara  = np.polyfit(self.PointTemp, self.PointEntropy, 11)
        self.fEntropy = np.poly1d(self.fEntropyPara)
        self.FitEntropy = self.fEntropy(self.FitTempRange)
        self.fEntropyR2  = r2(self.PointEntropy, self.fEntropy(self.PointTemp))
        
        self.fEnthalpyPara = np.polyfit(self.PointTemp, self.PointEnthalpy, 5)
        self.fEnthalpy = np.poly1d(self.fEnthalpyPara)
        self.FitEnthalpy = self.fEnthalpy(self.FitTempRange)
        self.fEnthalpyR2  = r2(self.PointEnthalpy, self.fEnthalpy(self.PointTemp))
        
        
    def plot_datas(self):
        import matplotlib.pyplot as plt        
        
        plt.plot(self.PointTemp, - self.PointTemp * self.PointEntropy, 'o',
                 label = '-Temp x Entropy, $-T*S(T)$ (eV/atom)')
        
        plt.plot(self.PointTemp,self.PointEnthalpy, 'o',label = 'Enthalpy, $H(T)$ (eV/atom)')
        plt.plot(self.FitTempRange, - self.FitEntropy * self.FitTempRange, '--', 
                 label = 'Fitted $-T*S(T)$ (eV) R2=%4.4f' % self.fEntropyR2)
        
        plt.plot(self.FitTempRange, self.FitEnthalpy, '--',
                 label = 'Fitted $H(T)$ (eV) R2=%4.4f' % self.fEnthalpyR2)
        plt.legend()
        plt.xlim([min(self.PointTemp), self.MaxTemp])
#         plt.ylim([-3, 1])#[ymin,ymax])
        plt.xlabel('Temperature (K)')
        plt.ylabel('Energy (eV/atom)')
