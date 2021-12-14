#!/bin/python 
##/usr/bin/env python
##
# coding: utf-8


## Functions for Amorphous Generator
import os,sys,random,datetime
# import commands
import JH_lib as jh
import numpy as np
import pandas as pd
from subprocess import check_output

##########################################
# Define distance calculation function
##########################################

def point_in_box_vector(point, box):
    ini2point=[point[j]-box['BoxOrigin'][j] for j in range(3)]
    if np.dot(ini2point, box['BoxVeca']) > 0 and  np.dot(ini2point, box['BoxVecb']) > 0 and     np.dot(ini2point, box['BoxVecc']) > 0 and point[2] > box['BoxOrigin'][2] and    point[2] < box['BoxOrigin'][2] +box['BoxVeca'][2]+box['BoxVecb'][2]+box['BoxVecc'][2] :
        return True
    else:        
        return False

def point_in_box_simple(point, box, unitcell):
    if box['BoxOrigin'][0] < point[0] < unitcell[0][0] and box['BoxOrigin'][1] < point[1] < unitcell[1][1] and box['BoxOrigin'][2] < point[2] < unitcell[2][2]:
        return True
    else:
        return False
    
def distance(a,b) :
    return sum([(x-y)**2.0 for x,y in zip(a,b)])**0.5 ;


def add_atom_in_box(totatom, Box, InsertAtomDF,OutputPath, unitcell):
    exAtomPosition = []
    minDistance = InsertAtomDF['radious'].sum()/InsertAtomDF.shape[0] *2.5
    for j in range(InsertAtomDF.shape[0]):
        label = InsertAtomDF.element.iloc[j]
        N_atom = InsertAtomDF.N_atom.iloc[j]
        for natom in range(N_atom):
            exAtomPosition.append([label,0,0,0,'T','T','T'])
    exAtomPositionDF = pd.DataFrame(exAtomPosition, columns=['label','x', 'y', 'z', 'rx','ry','rz'])
    
    for newatom in range(InsertAtomDF.N_atom.sum()):
        tot_attempt = 0
        condition = True
        while  condition  and tot_attempt < InsertAtomDF.N_atom.sum()* 1000 : 
            NewAtomPosition0 = [ random.random() for i in range(3)]
            while min(np.dot (NewAtomPosition0, unitcell)) < minDistance/2 and  tot_attempt < InsertAtomDF.N_atom.sum()* 1000:
                NewAtomPosition0 = [ random.random() for i in range(3)]
            NewAtomPosition = np.dot (NewAtomPosition0, unitcell)
            tot_attempt = tot_attempt + 1
            exAtomPositionDF['distance'] =( (exAtomPositionDF['x'] - NewAtomPosition[0])**2+
                                            (exAtomPositionDF['y'] - NewAtomPosition[1])**2+
                                            (exAtomPositionDF['z'] - NewAtomPosition[2])**2 )**0.5
            condition = exAtomPositionDF['distance'].min() < minDistance
        
#         print(newatom, NewAtomPosition, exAtomPositionDF['distance'].min() >  minDistance, min(NewAtomPosition))
        exAtomPositionDF['x'].iloc[newatom] = NewAtomPosition[0]
        exAtomPositionDF['y'].iloc[newatom] = NewAtomPosition[1]
        exAtomPositionDF['z'].iloc[newatom] = NewAtomPosition[2]
    NewPositions = exAtomPositionDF.iloc[:,0:7].values.tolist()
    NewCompound = jh.component_from_positions(NewPositions)
    jh.w_poscar(NewPositions, compound = NewCompound, filename = OutputPath, 
                unitcell = unitcell, Selective = True)
    return True

def build_SprayinBox( OutputPath='./outPOSCAR.vasp', 
                           AtomDensity = 2.65,  N_Atoms = [1 ,2] , AtomName = ['Si', 'O'] , MaxAtom = 200,
                           potdir = '/vasp/POTCAR/PAW_PBE', InitCell=[[],[],[]]):
    ##########################################
    # Preallocate & grep MASS/Radious from POTCAR
    ##########################################
    atommass = []
    atomradious = []
    for atom in AtomName:
        #Atomic Mass 
        atommass.append(atom_mass_radi_lib(atom, key='mass'))
        #Wigner-Seitz Radious 
        atomradious.append(atom_mass_radi_lib(atom, key='radi'))
    N_Atoms = np.array(N_Atoms) * (int(MaxAtom/sum(N_Atoms)))
    InsertAtomDF= pd.DataFrame(np.array([AtomName, N_Atoms,atommass,atomradious]).T,
                               columns=['element', 'N_atom', 'mass', 'radious'])

    InsertAtomDF['mass'] = InsertAtomDF['mass'].astype('float')
    InsertAtomDF['radious'] = InsertAtomDF['radious'].astype('float')
    InsertAtomDF['N_atom'] = InsertAtomDF['N_atom'].astype('int')#, 'N_atom':int, 'mass':float, 'radious': float) 
    InsertAtomDF['sumMass'] = InsertAtomDF['N_atom'] * InsertAtomDF['mass'] 


    ##########################################
    # Preallocate & for atom positions
    ##########################################
    totatom = sum(N_Atoms)
    ##########################################
    # Define Box for Inserting atom
    ##########################################
    Box = box_generator(InsertAtomDF,TargetDensity=AtomDensity, a=InitCell[0], b=InitCell[1], c=InitCell[2])  
    unitcell= [Box[ 'BoxVeca'], Box[ 'BoxVecb'] , Box[ 'BoxVecc']] 
#     add_atom_in_box(totatom,cellpar, Box, InsertAtomDF,OutputPath, unitcell)
    add_atom_in_box(totatom, Box, InsertAtomDF,OutputPath, unitcell)
    return True


def box_generator(InsertAtomDF,TargetDensity, a=[], b=[], c=[]):
    NumberofAvogadro = 6.022e23 # atom/mol
    Mass     = InsertAtomDF['sumMass'].sum()/NumberofAvogadro #g.atom/mol / (atom/mol ) = g
    Volume = Mass / TargetDensity * 1e24  # g/(g/cm^3) = cm^3 * *1e8)^3 = A^3
    
    if len(a)+len(b)+len(c) == 0:
        Height  = Volume**(1/3)
        Box ={'BoxOrigin':[0,0,0],    'BoxVeca': [Height,0,0],  'BoxVecb': [0,Height,0], 'BoxVecc':  [0,0,Height]}
    else:
        area = np.linalg.norm(np.cross(a,b))
        Height = Volume/area
        Box ={'BoxOrigin':[0,0,0],    'BoxVeca': a,  'BoxVecb': b, 'BoxVecc':  [0,0,Height]}
    return Box

def atom_mass_radi_lib(atom, key='radi'):
    #Lib: AtomicNumber,Label,Name,Radious,Mass
    lib= [[1,'H','Hydrogen',53,1.00], [2,'He','Helium',31,4.00], [3,'Li','Lithium',167,6.94], [4,'Be','Beryllium',112,9.01], [5,'B','Boron',87,10.81], 
    [6,'C','Carbon',67,12.01], [7,'N','Nitrogen',56,14.00], [8,'O','Oxygen',48,15.99], [9,'F','Fluorine',42,18.99], [10,'Ne','Neon',38,20.17],
    [11,'Na','Sodium',190,22.98], [12,'Mg','Magnesium',145,24.30], [13,'Al','Aluminium',118,26.98], [14,'Si','Silicon',111,28.08], [15,'P','Phosphorus',98,30.97], 
    [16,'S','Sulfur',88,32.06], [17,'Cl','Chlorine',79,35.45], [18,'Ar','Argon',71,39.09], [19,'K','Potassium',243,39.94], [20,'Ca','Calcium',194,40.08],
    [21,'Sc','Scandium',184,44.95], [22,'Ti','Titanium',176,47.90], [23,'V','Vanadium',171,50.94], [24,'Cr','Chromium',166,51.99], [25,'Mn','Manganese',161,54.93],
    [26,'Fe','Iron',156,55.84], [27,'Co','Cobalt',152,58.70], [28,'Ni','Nickel',149,58.93], [29,'Cu','Copper',145,63.54], [30,'Zn','Zinc',142,65.38], 
    [31,'Ga','Gallium',136,69.72], [32,'Ge','Germanium',125,72.59], [33,'As','Arsenic',114,74.92], [34,'Se','Selenium',103,78.96], [35,'Br','Bromine',94,79.90], 
    [36,'Kr','Krypton',88,83.80], [37,'Rb','Rubidium',265,85.46], [38,'Sr','Strontium',219,87.62], [39,'Y','Yttrium',212,88.90], [40,'Zr','Zirconium',206,91.22], 
    [41,'Nb','Niobium',198,92.90], [42,'Mo','Molybdenum',190,95.94], [43,'Tc','Technetium',183,98.00], [44,'Ru','Ruthenium',178,101.07], [45,'Rh','Rhodium',173,102.90],
    [46,'Pd','Palladium',169,106.40], [47,'Ag','Silver',165,107.86], [48,'Cd','Cadmium',161,112.41], [49,'In','Indium',156,114.82], [50,'Sn','Tin',145,118.69], 
    [51,'Sb','Antimony',133,121.75], [52,'Te','Tellurium',123,126.90], [53,'I','Iodine',115,127.60], [54,'Xe','Xenon',108,131.30], [55,'Cs','Cesium',298,132.90],
    [56,'Ba','Barium',253,137.33], [57,'La','Lanthanum',195,138.90], [58,'Ce','Cerium',185,140.12], [59,'Pr','Praseodymium',247,140.90], [60,'Nd','Neodymium',206,144.24], 
    [61,'Pm','Promethium',205,145.00], [62,'Sm','Samarium',238,150.40], [63,'Eu','Europium',231,151.96], [64,'Gd','Gadolinium',233,157.25], [65,'Tb','Terbium',225,158.92], 
    [66,'Dy','Dysprosium',228,162.50], [67,'Ho','Holmium',226,164.93], [68,'Er','Erbium',226,167.26], [69,'Tm','Thulium',222,168.93], [70,'Yb','Ytterbium',222,173.04], 
    [71,'Lu','Lutetium',217,174.96], [72,'Hf','Hafnium',208,178.49], [73,'Ta','Tantalum',200,180.94], [74,'W','Tungsten',193,183.85], [75,'Re','Rhenium',188,186.20],
    [76,'Os','Osmium',185,190.20], [77,'Ir','Iridium',180,192.22], [78,'Pt','Platinum',177,195.09], [79,'Au','Gold',174,196.96], [80,'Hg','Mercury',171,200.59], 
    [81,'Tl','Thallium',156,204.37], [82,'Pb','Lead',154,207.20], [83,'Bi','Bismuth',143,208.98], [84,'Po','Polonium',135,209.00], [85,'At','Astatine',127,210.00], 
    [86,'Rn','Radon',120,222.00], [87,'Fr','Francium','None',223.00], [88,'Ra','Radium','None',226.02], [89,'Ac','Actinium',195,227.02], [90,'Th','Thorium',180,231.03], 
    [91,'Pa','Protactinium',180,232.03], [92,'U','Uranium',175,237.04], [93,'Np','Neptunium',175,238.02], [94,'Pu','Plutonium',175,242.00], [95,'Am','Americium',175,243.00] ] 
    dicList=[]
    for i in lib:
        dicList.append({'label':i[1], 'radi':i[0], 'fullname':i[2], 'radious':i[3], 'mass':i[4]})
    elementDF = pd.DataFrame(dicList)
    if key == 'radi':
        return elementDF[elementDF.label == atom].radious.iloc[0] /100 # pm to Angtrom
    elif key == 'mass':
        return elementDF[elementDF.label == atom].mass.iloc[0]
    
from  tkinter import  Button, Label, Tk,  Grid, Entry, filedialog, END, StringVar,Text
from os import  listdir
## Functions for GUI


def readAllEnrty():
    workspace = OutputEntry.get()
    ElementList = [i.replace(' ','') for i in ElementEntry.get().split(',')]
    NumElementList =  [int(i.replace(' ','')) for i in NumElementEntry.get().split(',')]
    Density = float(DensityEntry.get())
    MaxAtomNum = int(MaxAtomNumEntry.get())
    NImage = int(NImageEntry.get())
    InitCell = [[float(UnitCellAxEntry.get()),float(UnitCellAyEntry.get()),float(UnitCellAzEntry.get())],
                [float(UnitCellBxEntry.get()),float(UnitCellByEntry.get()),float(UnitCellBzEntry.get())],
                [float(UnitCellCxEntry.get()),float(UnitCellCyEntry.get()),float(UnitCellCzEntry.get())]]

    outputText.delete('1.0', END)
    outputText.insert(END, str('Output Path/File: %s/POSCAR_XX\n'%workspace))
    
    
    if len(ElementList) == len(NumElementList):
        compisition=''
        Ndivision =0
        for element, Nelement in zip(ElementList, NumElementList):
            compisition+='%s%i' %(element,Nelement)
            Ndivision +=  Nelement

        outputText.insert(END, str('Compisition: %s\n'%compisition))
        outputText.insert(END, str('Density: %3.2f\n'%Density))
        outputText.insert(END, str('Maximun number of Atom: %i\n'%(int(MaxAtomNum/Ndivision) *int(Ndivision))))
        outputText.insert(END, str('Input Cell : %s\n'%(InitCell)))
        outputText.insert(END, str('Number of amorphous structure to build: %i\n'%(NImage)))
        
        outputText.insert(END, str('Start to construct the structure...\n\n: '))
        outputText.insert(END, str('In the Outputpath we generated\n: '))
        for i in range(NImage ):
            buiding = build_SprayinBox( OutputPath=workspace+'/POSCAR%2.2i.vasp' %i,
                                       N_Atoms = NumElementList , AtomName = ElementList, 
                                       MaxAtom = int(MaxAtomNum/Ndivision) *int(Ndivision),
                                   AtomDensity = Density,  
                            InitCell = InitCell
                                      )
            for n, i in enumerate(listdir(workspace)):
                if n % 3 == 0 and n !=0 :
                        outputText.insert(END, str('%s\n ' %i))
                else:   outputText.insert(END, str('%s  , ' %i))
    else:
        outputText.insert(END, str('Number of Elements: %s\n'%ElementList))
        outputText.insert(END, str('Wrong Input, Match the number of elements and component\n'))


        
def output(Inentry):
    path = filedialog.askdirectory()
    Inentry.delete(1, END)  # Remove current text in entry
    Inentry.insert(0, path)  # Insert the 'path'
    return path


class MainApplication():
    def __init__(self, master):
        self.master = master
        self.master.title("Amorphous Builder")
#         label = Label(self.master, text="Test Callback", )

    def LabelEntry(self, StringValue, ColNum=0, RowNum=0, initialString=False):
        def _on_click(event):
            event.widget.delete(0, END)
        label = Label(self.master, text=StringValue).grid( column = ColNum, row = RowNum, pady=5, padx=5)
        entry = Entry(self.master, width=40)
#         if initialString: 
        entry.grid( column = ColNum+1, row = RowNum, sticky='W', pady=5, padx=5)
        entry.bind("<Button-1>", _on_click)
        return label, entry
    
    
    def BrowsDirButton(self, entry, ColNum=0, RowNum=0):
        button = Button(self.master, text="Browse", command=lambda: output(entry))
        button.grid(column = ColNum, row = RowNum) 
        return button

    def InitCellInforTable(self, StringValue, ColNum=0, RowNum=0, initialString=False):
        def _on_click(event):
            event.widget.delete(0, END)
        axlabel = Label(self.master, text=StringValue).grid( column = ColNum, row = RowNum, pady=5, padx=5)
        xentry = Entry(self.master, width=10)#         if initialString: 
        yentry = Entry(self.master, width=10)#         if initialString: 
        zentry = Entry(self.master, width=10)#         if initialString: 
        xentry.grid( column = ColNum+1, row = RowNum, sticky='W', pady=5, padx=5)
        yentry.grid( column = ColNum+1, row = RowNum, sticky='W', pady=5, padx=5+30*3)
        zentry.grid( column = ColNum+1, row = RowNum, sticky='W', pady=5, padx=5+30*6)

        return axlabel,xentry,yentry,zentry

    def close(self):
        self.master.quit()
        return
    
    
    
if __name__ == '__main__':
    root = Tk()
    gui = MainApplication(root)
    rown = 0 
    IntroductionLabel=Label(gui.master, text='(POSCAR Format / based on Spraying)', font="Helvetica 12 bold").grid(column = 1, row = rown, columnspan =2 , pady=5, padx=5)
    
    rown +=1    
    InLabel = Label(gui.master, text='Type Input values', font="Helvetica 12 bold").grid(column = 1, row = rown, columnspan =2 , pady=5, padx=5)
    rown +=1
    ElementLabel,ElementEntry = gui.LabelEntry('Element', 1, rown)
    ElementExample=StringVar(gui.master, value='(example for Si3N4) Si, N')
    ElementEntry.configure(textvariable=ElementExample); rown +=1

    NumElementLabel,NumElementEntry = gui.LabelEntry('Number of Element', 1, rown)
    NumElementExample=StringVar(gui.master, value='(example for Si3N4) 3, 4')
    NumElementEntry.configure(textvariable=NumElementExample); rown +=1

    DensityLabel,DensityEntry = gui.LabelEntry('Density (g/cm^3)', 1, rown); rown +=1

    NImageLabel,NImageEntry = gui.LabelEntry('Number of Structure', 1, rown); rown +=1

    MaxAtomNumLabel,MaxAtomNumEntry = gui.LabelEntry('Max Atom Number', 1, rown); rown +=1

    OutputLabel,OutputEntry = gui.LabelEntry('Output Path', 1, rown)
    OuyputBrowse2 = gui.BrowsDirButton(OutputEntry, 3, rown); rown +=1

    defLabel = Label(gui.master, text='Initial Cell Condition', font="Helvetica 12 bold")
    defLabel.grid(column = 1, row = rown, columnspan =2 , pady=5, padx=5); rown +=1
    
    UnitCellALabel,UnitCellAxEntry,UnitCellAyEntry,UnitCellAzEntry = gui.InitCellInforTable('a-axis', 1, rown); rown +=1
    UnitCellBLabel,UnitCellBxEntry,UnitCellByEntry,UnitCellBzEntry = gui.InitCellInforTable('b-axis', 1, rown); rown +=1
    UnitCellCLabel,UnitCellCxEntry,UnitCellCyEntry,UnitCellCzEntry = gui.InitCellInforTable('c-axis', 1, rown); rown +=1
    
    UnitCellCxEntryExp=StringVar(gui.master, value='0')
    UnitCellCxEntry.configure(textvariable=UnitCellCxEntryExp)        
    UnitCellCyEntryExp=StringVar(gui.master, value='0')
    UnitCellCyEntry.configure(textvariable=UnitCellCyEntryExp)
    UnitCellCzEntryExp=StringVar(gui.master, value='0')
    UnitCellCzEntry.configure(textvariable=UnitCellCzEntryExp)
    
    defLabel2 = Label(gui.master, text='(put 0,0,0 in c-axis) for automatic vertical cell define')
    defLabel2.grid(column = 1, row = rown, columnspan =2 , pady=5, padx=5); rown +=1

    begin_button = Button(gui.master, text='Begin!', command=lambda: readAllEnrty())
    begin_button.grid(column = 1, row= rown+1, columnspan=3); rown +=2


    outputLabel = Label(gui.master, text='Progress Box', font="Helvetica 10 bold")
    outputLabel.grid(column = 1, row = rown, columnspan =1 , pady=5, padx=5,  sticky='w' ); rown +=1
    outputText=Text(gui.master,height=10)#,warp='word')
    outputText.grid(column = 1, row= rown, columnspan=3)

    root.mainloop()

    
