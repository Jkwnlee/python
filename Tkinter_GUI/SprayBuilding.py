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
#         print(True, point, box['BoxOrigin'] )
        return True
    else:        
#         print(False, point, box['BoxOrigin'] )
        return False

def point_in_box_simple(point, box, unitcell):
    if box['BoxOrigin'][0] < point[0] < unitcell[0][0] and        box['BoxOrigin'][1] < point[1] < unitcell[1][1]and        box['BoxOrigin'][2] < point[2] < unitcell[2][2]:
        return True
    else:
        return False
    
def distance(a,b) :
    return sum([(x-y)**2.0 for x,y in zip(a,b)])**0.5 ;

def add_atom_in_box(totatom,cellpar, Box,InsertAtomDF,OutputPath, unitcell):
    exAtomPosition = [[0 for i in range(3)] for j in range(totatom)]
    N_Atoms = InsertAtomDF['N_atom'].values.tolist()
    AtomName=InsertAtomDF['element'].values.tolist()
    AtomRadious=InsertAtomDF['radious'].values.tolist()
    LabelN = 0; maxAtom = N_Atoms[LabelN]
    tot_attempt = 0;
    NewPositions =[]
    for newatom in range(totatom):
        exatom = -1;
        while exatom < newatom and tot_attempt < totatom*100:
            NewAtomPosition = [ cellpar[i] * random.random() for i in range(3)] 
            tot_attempt = tot_attempt + 1
            
            if tot_attempt > totatom*100:     # Exit Loop if it takes too long
                print ("[CODE] WARNING!!! The LOOP TAKE TOO LONG"); break
                
            if point_in_box_vector(NewAtomPosition,Box ) == True :
                i = 0
                x = exAtomPosition[i]
                while distance(NewAtomPosition,x) > min(AtomRadious)*2.5 and i < totatom:
                    i = i + 1
                    if exAtomPosition[i] !=[0, 0, 0]: x = exAtomPosition[i]
                    else: break
                if distance(NewAtomPosition,x) > min(AtomRadious)*2.5: exatom = totatom
                else: pass
                
                if distance(NewAtomPosition,x) < 3:
                    pass
#                     print(distance(NewAtomPosition,x))
        
        if newatom < maxAtom:pass
        elif newatom ==  maxAtom:
            LabelN= LabelN + 1
            maxAtom = N_Atoms[LabelN] + maxAtom
        NewPositions.append([AtomName[LabelN],
                             NewAtomPosition[0],NewAtomPosition[1],NewAtomPosition[2], 
                             'T','T','T' ])
        
    NewCompound = jh.component_from_positions(NewPositions)
    jh.w_poscar(NewPositions, compound = NewCompound, filename = OutputPath, 
                unitcell = unitcell, Selective = True)
    return True

def build_SprayinBox( OutputPath='./outPOSCAR.vasp', 
                           AtomDensity = 2.65,  N_Atoms = [1 ,2] , AtomName = ['Si', 'O'] , MaxAtom = 200,
                           potdir = '/vasp/POTCAR/PAW_PBE'):
    ##########################################
    # Preallocate & grep MASS/Radious from POTCAR
    ##########################################
    atommass = []
    atomradious = []
    for atom in AtomName:
        #Atomic Mass 
        batcmd = ''.join(['grep MASS ',potdir,'/',atom,'/POTCAR | awk \'{print $3}\' | cut -d";" -f1'])
        atommass.append(float(check_output(batcmd, shell = True)))
        #Wigner-Seitz Radious 
        batcmd = ''.join(['grep RWIGS ',potdir,'/',atom,'/POTCAR | awk \'{print $3}\' | cut -d";" -f1'])
        atomradious.append(float(check_output(batcmd, shell = True)))
    N_Atoms = np.array(N_Atoms) * (int(MaxAtom/sum(N_Atoms)))
#     print(N_Atoms,MaxAtom, int(MaxAtom/sum(N_Atoms)) , N_Atoms, list(N_Atoms))
    
    InsertAtomDF= pd.DataFrame(np.array([AtomName, N_Atoms,atommass,atomradious]).T,
                               columns=['element', 'N_atom', 'mass', 'radious'])

    InsertAtomDF['mass'] = InsertAtomDF['mass'].astype('float')
    InsertAtomDF['radious'] = InsertAtomDF['radious'].astype('float')
    InsertAtomDF['N_atom'] = InsertAtomDF['N_atom'].astype('int')#, 'N_atom':int, 'mass':float, 'radious': float) 
    InsertAtomDF['sumMass'] = InsertAtomDF['N_atom'] * InsertAtomDF['mass'] 

    NumberofAvogadro = 6.022e23 # atom/mol
    Mass     = InsertAtomDF['sumMass'].sum()/NumberofAvogadro #g.atom/mol / (atom/mol ) = g
    Volume = Mass / AtomDensity * 1e24  # g/(g/cm^3) = cm^3 * *1e8)^3 = A^3
    Height  = Volume**(1/3)
    
    ##########################################
    # Preallocate & for atom positions
    ##########################################
    totatom = sum(N_Atoms)
    ##########################################
    # Define Box for Inserting atom
    ##########################################
    Box ={'BoxOrigin':[0,0,0],    'BoxVeca': [Height,0,0],  'BoxVecb': [0,Height,0], 'BoxVecc':  [0,0,Height],}
    cellpar = [Height,Height,Height]
    unitcell= [Box[ 'BoxVeca'], Box[ 'BoxVecb'] , Box[ 'BoxVecc']]
    add_atom_in_box(totatom,cellpar, Box, InsertAtomDF,OutputPath, unitcell)
    return True



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
    PotcarPathIn=PotcarPathEntry.get()


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
        outputText.insert(END, str('Number of amorphous structure to build: %i\n'%(NImage)))

        for i in range(NImage ):
            buiding = build_SprayinBox( OutputPath=workspace+'/POSCAR%2.2i' %i,
                                       N_Atoms = NumElementList , AtomName = ElementList, 
                                       MaxAtom = int(MaxAtomNum/Ndivision) *int(Ndivision),
                                   AtomDensity = Density,  
                                   potdir = PotcarPathIn)
        outputText.insert(END, str('In the Outputpath we generated\n: '))
        for n, i in enumerate(listdir(workspace)):
            if n % 5 == 0 and n !=0 :   outputText.insert(END, str('%s\n ' %i))
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
    ElementExample=StringVar(gui.master, value='example for Si3N4) Si, N')
    ElementEntry.configure(textvariable=ElementExample); rown +=1

    NumElementLabel,NumElementEntry = gui.LabelEntry('Number of Element', 1, rown)
    NumElementExample=StringVar(gui.master, value='example for Si3N4) 3, 4')
    NumElementEntry.configure(textvariable=NumElementExample); rown +=1

    DensityLabel,DensityEntry = gui.LabelEntry('Density (g/cm^3)', 1, rown); rown +=1

    NImageLabel,NImageEntry = gui.LabelEntry('Number of Structure', 1, rown); rown +=1

    MaxAtomNumLabel,MaxAtomNumEntry = gui.LabelEntry('Max Atom Number', 1, rown); rown +=1

    OutputLabel,OutputEntry = gui.LabelEntry('Output Path', 1, rown)
    OuyputBrowse2 = gui.BrowsDirButton(OutputEntry, 3, rown); rown +=1

    defLabel = Label(gui.master, text='Default is set below', font="Helvetica 12 bold")
    defLabel.grid(column = 1, row = rown, columnspan =2 , pady=5, padx=5); rown +=1

    PotcarPathLabel,PotcarPathEntry = gui.LabelEntry('POTCAR Path (for Mass)', 1, rown)
    PotcarPath=StringVar(gui.master, value='~/vasp/POTCAR/PAW_PBE')
    PotcarPathEntry.configure(textvariable=PotcarPath)
    PotcarBrowse2 = gui.BrowsDirButton(PotcarPathEntry, 3, rown); rown +=1

    begin_button = Button(gui.master, text='Begin!', command=lambda: readAllEnrty())
    begin_button.grid(column = 1, row= rown+1, columnspan=3); rown +=2


    outputLabel = Label(gui.master, text='Progress Box', font="Helvetica 10 bold")
    outputLabel.grid(column = 1, row = rown, columnspan =1 , pady=5, padx=5,  sticky='w' ); rown +=1
    outputText=Text(gui.master,height=10)#,warp='word')
    outputText.grid(column = 1, row= rown, columnspan=3)

    root.mainloop()





