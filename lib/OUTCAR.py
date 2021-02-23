
class OUTCAR():
    outlist = []
    SysName =''
    Spin= 0
    TargetText = ''
    outcartexts = []
    fermilevels = []
    energylist  = []
    
    #merge from DOS class
    NKpts = 0	# Number of K-points
    NBands = 0	# Number of bands
    
    DGap = 0	# Smallest direct gap
    DGapKp = 0	# Kpoint with smalles direct gap
    VBMStKp = 0	# Kpoint with highest occupied state
    CBMStKp = 0	# Kpoint with lowest unoccupied state
    EVBM = 0	# Highest occupied state in Global 
    ECBM = 0	# Lowest occupied state in Global
    EGap = 0		# Global gap (might be indirect)
    
    UP_DGap = 0	# Smallest direct gap
    UP_DGapKp = 0	# Kpoint with smalles direct gap
    UP_VBMStKp = 0	# Kpoint with highest occupied state
    UP_CBMStKp = 0	# Kpoint with lowest unoccupied state
    UP_EVBM = 0	# Highest occupied state in Global 
    UP_ECBM = 0	# Lowest occupied state in Global
    UP_EGap = 0		# Global gap (might be indirect)
    
    DOWN_DGap = 0	# Smallest direct gap
    DOWN_DGapKp = 0	# Kpoint with smalles direct gap
    DOWN_VBMStKp = 0	# Kpoint with highest occupied state
    DOWN_CBMStKp = 0	# Kpoint with lowest unoccupied state
    DOWN_EVBM = 0	# Highest occupied state in Global 
    DOWN_ECBM = 0	# Lowest occupied state in Global
    DOWN_EGap = 0		# Global gap (might be indirect)
    
    Met = False	# Switches to true if the spin channel is metallic (partial occupations)
    
    def __init__(self, filename):
        self.outcar = open(filename, 'r')
        self.outcartexts = open(filename, 'r').readlines()
        self.OutcarText = open(filename, 'r').read()
        
        #merge from DOS class
        self.LocVBSt=-10000000	# Setting low initial limit for highest occupied state	
        self.LocCBSt=10000000	# Setting high initial limit for lowest unoccupied state
        self.VBSt  = []		# Array for local highest occupied state
        self.CBSt  = []		# Array for local lowest unoccupied state
        self.D_Gaps = []	# Array for direct gaps
        
        self.UP_VBSt  = []		# Array for local highest occupied state
        self.UP_CBSt  = []		# Array for local lowest unoccupied state
        self.UP_D_Gaps = []	# Array for direct gaps
        self.DOWN_VBSt  = []		# Array for local highest occupied state
        self.DOWN_CBSt  = []		# Array for local lowest unoccupied state
        self.DOWN_D_Gaps = []	# Array for direct gaps
    def get_fermi(self):
        return [float(x.split()[2]) for x in self.outcartexts if 'E-fermi' in x][-1]
    
    def get_energy(self):
        return [float(x.split()[4]) for x in self.outcartexts if 'energy without entropy' in x][-1]
    
    def get_nelect(self):
        return int([float(x.split()[2]) for x in self.outcartexts if 'NELECT =' in x][-1])
    
    def get_system(self):
        return [ x.split()[2] for x in self.outcartexts if "SYSTEM =" in x][-1]
    
    def get_spin(self):
        return [ int(x.split()[2]) for x in self.outcartexts if "ISPIN  =" in x][-1]

    def get_1sCoreLevel(self):
        return [ x.split()[2] for x in self.outcartexts if "  1s  " in x]
       
    def get_LineNr(self):
        if 'k-point     1 :' in self.OutcarText: # VASP has switched Output format between 5.2 and 5.3,
            TargetText='k-point     1 :'	# this should now work in both (all?) cases
        else:
            TargetText="k-point   1"
        return [ i for i, x in enumerate(self.outcartexts) if TargetText in x]
   
    def get_NKPTS(self):
        outlist = [ x.split() for x in self.outcartexts if "NKPTS =" in x][-1]
        return int(outlist[14]), int(outlist[3])
    
    def get_bandgap(self):
        NonMet=[1.0, 2.0, 0.0]	# All other occupation numbers indicate a metallic state
        Occ=[1.0,2.0]		# Possible occcupation numbers for full occupation
        Data=[] # Initializing array for bandenergy and occupation data
        Spin = self.get_spin()
        LineNr = self.get_LineNr()
        Nbands, Nkpts = self.get_NKPTS()
        if Spin == 1:
            StartLine=LineNr[-1]
            EndLine=StartLine + Nkpts * Nbands + Nkpts*3
        elif Spin == 2:
            StartLine=LineNr[-2]
            EndLine=StartLine + 2*Nkpts*Nbands+ 2*Nkpts*3
        else:
            print( " " )
            print( " Error encountered regarding spin states. Aborting. " )
            print( " " )
            
        not_these=["k-point", "band", "spin"]
        #filling data array
        for count, line in enumerate(self.outcartexts,1):
            if not line == '\n':	# Discarding empty lines
                if count >= StartLine and count <= EndLine and line.split()[0] not in not_these:
                    Data.append([float(line.split()[1]), float(line.split()[2])])

        if Spin  == 1:
            for k in range(0,Nkpts):
                for b in range(0,Nbands):
                    if Data[k*Nbands+b][1] not in NonMet:
                        self.Met=True
                    elif Data[k*Nbands+b][1] in Occ and Data[k*Nbands+b][0] > self.LocVBSt:
                        self.LocVBSt=Data[k*Nbands+b][0]
                    elif Data[k*Nbands+b][1] == 0.0 and Data[k*Nbands+b][0] < self.LocCBSt:
                        self.LocCBSt=Data[k*Nbands+b][0]
                self.D_Gaps.append(self.LocCBSt-self.LocVBSt)
                self.VBSt.append(self.LocVBSt)
                self.CBSt.append(self.LocCBSt)
                self.LocVBSt = -10000000
                self.LocCBSt =  10000000
            self.DGap, self.DGapKp = min((self.DGap, self.DGapKp) for (self.DGapKp, self.DGap) in enumerate(self.D_Gaps,1))
            self.CBStKp = self.CBSt.index(min(self.CBSt))+1
            self.VBStKp = self.VBSt.index(max(self.VBSt))+1
            self.EVBM= max(self.VBSt)
            self.ECBM= min(self.CBSt)
            self.EGap = self.ECBM-self.EVBM
            
        elif Spin  ==2:
            Up_Data=Data[:len(Data)/2]
            for k in range(0,Nkpts):
                for b in range(0,Nbands):
                    if Up_Data[k*Nbands+b][1] not in NonMet:
                        self.Met=True
                    elif Up_Data[k*Nbands+b][1] in Occ and Up_Data[k*Nbands+b][0] > self.LocVBSt:
                        self.LocVBSt=Data[k*Nbands+b][0]
                    elif Up_Data[k*Nbands+b][1] == 0.0 and Up_Data[k*Nbands+b][0] < self.LocCBSt:
                        self.LocCBSt=Data[k*Nbands+b][0]
                self.UP_VBSt.append(self.LocVBSt)
                self.UP_CBSt.append(self.LocCBSt)
                self.UP_D_Gaps.append(self.LocCBSt-self.LocVBSt)
                self.LocVBSt = -10000000
                self.LocCBSt =  10000000
            self.UP_DGap, self.UP_DGapKp = min((self.UP_DGap, self.UP_DGapKp) for (self.UP_DGapKp, self.UP_DGap) in enumerate(self.UP_D_Gaps,1))
            self.UP_CBStKp = self.UP_CBSt.index(min(self.UP_CBSt))+1
            self.UP_VBStKp = self.UP_VBSt.index(max(self.UP_VBSt))+1
            self.UP_EVBM= max(self.UP_VBSt)
            self.UP_ECBM= min(self.UP_CBSt)
            self.UP_EGap = self.UP_ECBM-self.UP_EVBM
            
            
            Down_Data=Data[len(Data)/2:]
            for k in range(0,Nkpts):
                for b in range(0,Nbands):
                    if Down_Data[k*Nbands+b][1] not in NonMet:
                        self.Met=True
                    elif Down_Data[k*Nbands+b][1] in Occ and Down_Data[k*Nbands+b][0] > self.LocVBSt:
                        self.LocVBSt=Data[k*Nbands+b][0]
                    elif Down_Data[k*Nbands+b][1] == 0.0 and Down_Data[k*Nbands+b][0] < self.LocCBSt:
                        self.LocCBSt=Data[k*Nbands+b][0]
                self.DOWN_VBSt.append(self.LocVBSt)
                self.DOWN_CBSt.append(self.LocCBSt)
                self.DOWN_D_Gaps.append(self.LocCBSt-self.LocVBSt)
                self.LocVBSt = -10000000
                self.LocCBSt =  10000000
            self.DOWN_DGap, self.DOWN_DGapKp = min((self.DOWN_DGap, self.DOWN_DGapKp) for (self.DOWN_DGapKp, self.DOWN_DGap) in enumerate(self.DOWN_D_Gaps,1))
            self.DOWN_CBStKp = self.DOWN_CBSt.index(min(self.DOWN_CBSt))+1
            self.DOWN_VBStKp = self.DOWN_VBSt.index(max(self.DOWN_VBSt))+1
            self.DOWN_EVBM= max(self.DOWN_VBSt)
            self.DOWN_ECBM= min(self.DOWN_CBSt)
            self.DOWN_EGap = self.DOWN_ECBM-self.DOWN_EVBM
