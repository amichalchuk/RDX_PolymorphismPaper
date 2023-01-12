# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 16:25:51 2022

@author: michaaal
"""


# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 19:30:29 2022

@author: michaaal
"""


# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 09:08:22 2020
@author: amichalc
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 21:55:40 2020
@author: amichalc
NOTES:
    ** At the minute where we take a negaitve of e.g. w-wt we sometimes get negative values where lower phonon is on the edge of the Gaussian and the target is e.g. the peak. We can now make this either = 0, or = *-1. The latter reflects up-pumping of the relative magnitude of scattering... 
    ** Need to think this thourhg if a better approach might be needed. E.g. turn DOS into list of [0,1] and up-pump that way instead. Then also do not need arbitrary lower cut-off.
outline:
    1. Normalize phonon bath to 6Z+Y
    2. Shock temperatures g(w,T) = g(w)*n(w,Phi_Sh)
    3. Project to doorways; g(w,T)/Z ==> g(w) =/= 0
    4.
"""
#
## A Michalchuk, C Morrison University of Edinburgh
##
##TO CHECK
#
#  Tolerance function at line 158. When -0.01 sometimes throws up error.. why?
#  Add intrinsic calculation of shock based on adiabatic compression and Dulong-Petit heat capacity
#  Add heat capacity calcs in Einstein model from DOS 
#
###########################################3
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import sys
import csv
import pylab
from warnings import warn
import scipy.integrate as integrate
from os import getcwd
from tkinter.filedialog import askopenfilenames, asksaveasfilename
from tkinter import Tk


class DensityofStates:
    def __init__(self,freq,gox):
        self.frequency = freq
        self.gox = []# gox
        self.gox_nonnorm = gox *1
        self.gox_orig = []
        #self.gox = gox#/np.sum(self.gox)        
        self.otone_freq = []
        self.otone_gox = [] #GoX defined as gox/otone_num                              
        self.otone_temp_gox = np.zeros(np.size(freq)) #This is the temperature generated gox

        self.proj_otone = np.zeros(np.size(freq))
        self.proj_otonein = np.zeros(np.size(freq))
        self.proj_twoPDoS = np.zeros(np.size(freq))
      
        self.twoPDoS = np.zeros(np.size(freq))  #this is the wavenumber-resolved tpdos
        self.tpd_temp = [] #this is the counter for tpdos accumulation
        
        self.ShockFreq = []
        self.Shockgox =[]
        self.ShockPopul =[]
        self.ShockAll = np.zeros(np.size(freq))
        
        self.VibEnergy = np.zeros(np.size(freq))
        
    def NormaliseDOS(self,WMax):
        NormNum=np.int(input('Normalise to [1] 3N; [2] 1; [3] Phonon Bath Vibs; [5] 6Z+Y/Vc '))
        NMol=np.int(input('Number of Molecules'))
        
        if NormNum==1:
            self.gox = self.gox_nonnorm
            self.gox_orig = self.gox_nonnorm
            intcut=0.01
            print('normalising the gox to 3N -- done')
            return intcut,NMol
        
        if NormNum==2:
            self.gox = self.gox_nonnorm / np.trapz(self.gox_nonnorm)
            self.gox_orig = self.gox_nonnorm / np.trapz(self.gox_nonnorm)
            intcut=0.0001
            print('normalising the gox to 1 -- done')
            return intcut,NMol
        
        if NormNum==3: #normalise to Nband < OMax and to corresponding calibration for > Omax
            PBVal = np.int(input('Number of phonon modes in bath? :'))
            WMaxID = np.where(np.isclose(self.frequency,WMax,atol=0.5))
            WMaxID = WMaxID[0][0]
            Normval = PBVal / integrate.simps(self.gox_nonnorm[0:WMaxID+1])
            self.gox = self.gox_nonnorm*Normval
            self.gox_orig = self.gox_nonnorm*Normval
            
            intcut=0.01
            print('normalising the gox to Phonon Bath -- done')
            return intcut, NMol #np.floor(PBVal/6)
        
        if NormNum==4: #normalise by Nband phonon bath < Omax and atom number  > Omax
            PBVal = np.int(input('Number of phonon modes in bath? :'))
            AtNum = np.int(input('Atom Number: ')) * 3
            WMaxID = np.where(np.isclose(self.frequency,WMax,atol=0.5))
            WMaxID = WMaxID[0][0]
            Normval1 = PBVal / integrate.simps(self.gox_nonnorm[0:WMaxID+1])
            self.gox = self.gox_nonnorm*Normval1
            self.gox_orig = self.gox_nonnorm*Normval1
            Normval2 = AtNum / integrate.simps(self.gox_nonnorm[WMaxID:])
            self.gox = self.gox_nonnorm*Normval1
            self.gox[WMaxID:] = self.gox_nonnorm[WMaxID:]*Normval2
            self.gox_orig = self.gox*1
            
            intcut=0.01
            print('normalising the gox to Phonon Bath -- done')
            return intcut, NMol #np.floor(PBVal/6)   
 
        if NormNum==5: #normalise to Nband/Vol < OMax as Fried and Ruggiero
            PBVal = np.int(input('Number of phonon modes in bath? 6(Z+Y) :'))
            VOL = np.float(input('UC Volume / A^3 :   '))
            WMaxID = np.where(np.isclose(self.frequency,WMax,atol=0.5))
            WMaxID = WMaxID[0][0]
            Normval = (PBVal/VOL) / integrate.simps(self.gox_nonnorm[0:WMaxID+1])
            self.gox = self.gox_nonnorm*Normval
            self.gox_orig = self.gox_nonnorm*Normval
            
            intcut=0.00001
            print('normalising the gox to Phonon Bath -- done')
            return intcut, NMol #np.floor(PBVal/6)
        
        
        
    def Overtone_DoS(self, otone_num,ShOn): 
        n=0
        self.otone_freq = self.frequency * otone_num   
        if ShOn==1:
            self.otone_gox = self.ShockAll/otone_num  
        else:
            self.otone_gox = self.gox/otone_num  
            
            
def BoseEinstein(Temp,freq):
    freqhz = freq*29979245800 #convert wavenumber to Hz
    denom = np.exp((hbar*freqhz)/(boltz*Temp))-1
    if denom <= 0:
        popul=0
        warn('Frequency < 0')
    else: 
        popul = 1/denom
    #print('Population at ', freqhz/29979245800, ' is ', popul)
    return popul 


def Temp(DoS,WMax,toler,ShOn):      ##This bit of code checked and seems OK. I don't like that it overwrites DoS.gox... but c'est la vie for l'instant 
   WmaxID=np.where(np.isclose(DoS.frequency,WMax,atol=toler))
   WmaxID=WmaxID[0][-1]
   #WmaxID=WmaxID[0][0]
   n=0
   if ShOn==1:
       ShTemp = float(input("Shock Temperatue: "))
   else:
       ShTemp=Eqm_Temp*1
  
   for obj in DoS.frequency:
       if obj <= WMax:
           DoS.ShockAll[n]=DoS.gox_orig[n]*BoseEinstein(ShTemp,obj)
           n+=1
       else:
           DoS.ShockAll[n]=DoS.gox_orig[n]*BoseEinstein(Eqm_Temp,obj)
           n+=1
   
   FigDOSALL=plt.figure('DOS Shock All')
   plt.plot(DoS.frequency, DoS.ShockAll, color='k',label=('shock all Wmax=%.1f'%WMax))
   plt.xlabel('Wavenumber /cm$^{-1}$',size=18)
   plt.xlim(0,2000)
   plt.xticks(size=14)
   plt.ylabel('$g(\omega) $ x $ n(T,\omega)$  / states.a.u.$^{-1}$',size=18)
   plt.ylim(0,np.max(DoS.ShockAll)+10)
   plt.yticks(size=14)
   vline1=plt.axvline(185,linestyle='--',color='black')
   plt.legend()
   
   return ShTemp

def Proj_DoS(DoS,XX,otone_num,PBVal):
    n1=0  
    if XX==1:    #OTONE
        FREQPROJ = DoS.otone_freq.copy() 
        DOSPROJ = DoS.proj_otonein                                  
    if XX==2:   #Two PDOS
        FREQPROJ = DoS.frequency.copy()
        DoS.proj_twoPDoS = DoS.twoPDoS.copy() 
        DOSPROJ= DoS.proj_twoPDoS
        
    for obj in DoS.gox_orig:
        ProjOTID = np.where(np.isclose(FREQPROJ,DoS.frequency[n1],atol=toler*otone_num+0.01))
        if isinstance(ProjOTID, int) == False:
            try:
                ProjOTID = ProjOTID[0][-1]
            except IndexError:
                print("Seems no ProjOTID!")
                print(" ")
                print(n1,'',toler*otone_num+0.01,'',DoS.frequency[n1])
                sys.exit(1)
            
        if obj <= intcut:               ##must choose here a good cutoff for minimum projection onto fundamentals
            DOSPROJ[ProjOTID] = 0
           # print('removed point at ', ProjOTID)
        else: 
            DOSPROJ[ProjOTID] = DOSPROJ[ProjOTID] / PBVal            #This now normalises the projection onto the number of states it falls into.
             
        n1=n1+1

    FigOTONEALL=plt.figure('Projected Overtone')
    if XX==1:
        plt.plot(DoS.otone_freq, DoS.proj_otonein, linestyle='--',label=('ProjOT Wmax=%.1f'%WMax))
    if XX==2:
        plt.plot(DoS.frequency,DoS.proj_twoPDoS, linestyle='-',label=('Proj2PDOS Wmax=%.1f'%WMax))
    plt.xlabel('Wavenumber /cm$^{-1}$',size=18)
    plt.xlim(0,2000)
    plt.xticks(size=14)
    plt.ylabel('$g(\omega) $ x $ n(T,\omega)$  / states.a.u.$^{-1}$',size=18)
    plt.ylim(0,np.max(DoS.proj_twoPDoS)+10)
    plt.yticks(size=14)
    vline1=plt.axvline(WMax,linestyle='--',color='black')
    plt.legend()

    return FigOTONEALL,DOSPROJ,intcut

def Make2Step(DoS,ShOn):
    otone_num = int(input('Overtone Order: '))
    DoS.Overtone_DoS(otone_num,ShOn)
    
    if "y" in input('Project Overtone onto Fundamental DoS?  '):
        #if ShOn==1:
        DoS.proj_otonein = DoS.otone_gox.copy() # DoS.ShockAll.copy()                        

        OTONEFIG,DoS.proj_otone,intcut = Proj_DoS(DoS,1,otone_num,PBVal)
        print('Otone projected onto fundamentals - Done')    

    WmaxID=np.where(np.isclose(DoS.frequency,WMax,atol=toler))
    WmaxID=WmaxID[0][-1]
    Wmax2ID=np.where(np.isclose(DoS.frequency,2*WMax,atol=toler))
    Wmax2ID=Wmax2ID[0][-1]
    WmaxID=int(WmaxID)

    print('Wmax 1', WmaxID)
    print('Wmax 2',Wmax2ID)
    
    TempProjOT = np.zeros(np.size(DoS.proj_otone)*2)
    itval=0
   # itorig=0
    for it in TempProjOT[:-5]:
        if itval/2 % 1 == 0.0:
            TempProjOT[itval] = DoS.proj_otone[np.int(itval/2)]
            itval+=1
            #itorig+=1
        else: #this creates an interpolation between data points because the otone projection is doubled in the x-axis!
            useit=((itval/2))
            TempProjOT[itval] = (DoS.proj_otone[int(useit+0.5)] + DoS.proj_otone[int(useit-0.5)])/2
            itval+=1
            #itorig+=1
    
    for val in np.arange(0,WMax,step=1):
        val=int(val)
        DoS.gox[WmaxID+val] = TempProjOT[WmaxID+val] - DoS.ShockAll[WmaxID+val]  # + shockall will be n_OT + n_eqm. -shockall will be n_OT - n_eqm
        if DoS.gox[WmaxID+val] < 0:
            DoS.gox[WmaxID+val] = DoS.gox[WmaxID+val]*(-1)#0
    
    
    indnum=0
    for ind in DoS.gox:
        
        if ind < 0:
            warn('-----WARNING: Values of Projected Overtones are NEGATIVE! -----')
            print(DoS.frequency[indnum])
        indnum+=1
    
  #  for repop in np.arange(0,np.size(DoS.gox[0:WmaxID])):
  #      repop=int(repop)
  #      DoS.gox[repop]=DoS.gox_orig[repop]*BoseEinstein(ShTemp,DoS.frequency[repop])
        
    
    plt.figure('DoSafter Step 2')
    plt.plot(DoS.frequency,DoS.gox,label='Step 2 DoS Wmax=%.1f'%WMax)
    plt.legend()
    return otone_num, intcut

def RevBE(vib,occ):
    if occ <= 0:
        Temp = 0
    else: 
        freqhz = vib*29979245800 #convert wavenumber to Hz
        Temp = hbar*freqhz/(np.log(1+(1/occ))*boltz) #np.exp((hbar*freqhz)/(boltz*Temp))-1
    #print(occ)
    return Temp

def NtoEnergy(DoS):
    n=0
    for vib in DoS.proj_twoPDoS:
        DoS.VibEnergy[n] = RevBE(DoS.frequency[n],vib)
        n+=1
    return np.sum(DoS.VibEnergy)


## USER INPUTS VARIABLES
wd=getcwd()
main_window = Tk()
fname = askopenfilenames(initialdir=wd, filetypes=(('dos files', '*.dos'),('xy files', '*.xy'), ('dat files', '*.dat')))
main_window.destroy()
fname=str(fname)
#fname = input('Density of states file: ')
Eqm_Temp = float(input('Equilibrium Temperature /K: '))
WMax = float(input('Omega Max /cm^-1: '))

## DEFINE Fundamental Constants ##################################################################################
hbar  = 1.054571800*10**(-34) ## J.s
boltz = 1.38064852*10**(-23) ## J.K-1

## Setup initial variabls ########################################################################################  
datalist = np.loadtxt(fname[2:-3])
frequency,gox,useless = datalist.transpose()


try:
    CHECK
except NameError:
    CHECK=1
if CHECK==0:
    print('adding some checkpoing data to output -- done')

DoS=DensityofStates(frequency,gox)    ###Original DoS from QM / MM file in x=freq, y=g(w)
intcut,PBVal=DoS.NormaliseDOS(WMax)
toler = ((DoS.frequency[3] - DoS.frequency[2])/2)+0.01  ##Tolerance on values due to finite frequency spacing##

###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
#                                                    MAIN CODE STARTS HERE
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################



if "y" in input('Include shock temperature mechanism?  '):
    ShOn=1
    print ("  --Shock Populations Will be Set--  ")
else: 
    ShOn=0
    print ("  --Shock Populations Will NOT be Set--  ")    
    
ShTemp=Temp(DoS,WMax,toler,ShOn)  ##Here we generate DoS.ShockAll

if "y" in input('Use two-step VUP mechanishm?  '):
    otone_num, intcut = Make2Step(DoS,ShOn)
    print('Projected Otone onto gox for 2-step model ' )

for w in DoS.frequency:  #Define target frequency w
    DoS.tpd_temp = 0
    freqID = np.where(DoS.frequency == w)
    freqID = freqID[0][0]
    freqminID=np.where(DoS.frequency <= 1.0)
    try:
        freqminID=freqminID[0][-1]
    except IndexError:
        print('Seems no freqminID!')
        print(' ')
        print('frequency= ', w)
    #if DoS.frequency[freqID] >= 2*WMax:                                           # This ensures omega1 is a phonon or doorway
    #    freqID_2W = np.where(DoS.frequency <= 2*WMax)
    #    freqID_2W = freqID_2W[0][-1]
    #    New_Frequency_List = DoS.frequency[0:freqID_2W]*1
    #else:
     #   New_Frequency_List = DoS.frequency[0:freqID]*1
    New_Frequency_List = DoS.frequency[freqminID:freqID]*1
    #print('omega is: ', w)
    
    for omega1 in New_Frequency_List[:-1]:
        
        if omega1<= WMax:
            omega2=omega1*1 #If omega1 is a phonon mode, then only overtone scattering happens
        else:
            omega2 = w - omega1 # if omega1 is doorway mode, then omega2 must be w-w_1
     
        if w >=WMax and w < 2*WMax and omega2 <=WMax: #if target is doorway and w2 is phonon
            #already have that only overtone contributes to this region. # and omega1 <=WMax:
          
                w2IDint = np.where(np.isclose(DoS.frequency,omega2,atol=toler))
                
                if isinstance(w2IDint, int) == False:
                    w2IDint=w2IDint[0][0]  #This gives us the element position of omega2 in the frequency array
                    
                w1IDint=np.where(DoS.frequency==omega1)
                w1ID=w1IDint[0][0] 
                
                
                w2ID=w2IDint
                GoX_omega2 = DoS.gox_orig[w2ID]
                GoX_omega1 = DoS.gox_orig[w1ID]
               
                if omega1 == omega2 and np.isclose(omega1,w/2,atol=toler) and GoX_omega1 > intcut :
                #If w1,w2 are both phonons, and are the otone of target
                    num_doesexist =  (GoX_omega1*0.5*BoseEinstein(ShTemp,omega1) - DoS.gox_orig[freqID]*BoseEinstein(Eqm_Temp,w))


                #Next two lines include self scattering inside the doorway modes

                elif omega1 > WMax and omega1 < 2*WMax and omega2 < WMax and GoX_omega1 > intcut and GoX_omega2 > intcut :
                    num_doesexist = (GoX_omega1*BoseEinstein(Eqm_Temp,omega1))+(GoX_omega2*BoseEinstein(ShTemp,omega2))-(DoS.gox_orig[freqID]*BoseEinstein(Eqm_Temp,w))  #(GoX_omega1+GoX_omega2 - DoS.gox_orig[freqID])*(BoseEinstein(Eqm_Temp,omega1)+BoseEinstein(Eqm_Temp,omega2)-BoseEinstein(Eqm_Temp,w))

                else:
                    num_doesexist = 0# (GoX_omega1*GoX_omega2) * (BoseEinstein(ShTemp,omega1)*BoseEinstein(ShTemp,omega2)-BoseEinstein(Eqm_Temp,w))
                 
                if num_doesexist <0:
                    num_doesexist = 0  
            
         
        elif w >= 2*WMax and omega2 <= WMax and omega1 >= WMax and omega1 <=2*WMax: 
            #Now we say if target is above doorway region, omega1 must be doorway, and omega2 is a phonon
                
                w2IDint = np.where(np.isclose(DoS.frequency,omega2,atol=toler))
                
                if isinstance(w2IDint, int) == False:
                    w2IDint=w2IDint[0][0]  #This gives us the element position of omega2 in the frequency array
                    
                w1IDint=np.where(DoS.frequency==omega1)
                w1ID=w1IDint[0][0] 
                
                w2ID=w2IDint
                GoX_omega2 = DoS.gox_orig[w2ID]
                GoX_omega1 = DoS.gox_orig[w1ID]
               
                if GoX_omega2 > intcut and GoX_omega1 > intcut:  
                    #Should we still keep phonon bath hot, or back to Eqm?
                   num_doesexist = DoS.gox[w1ID]+GoX_omega2*BoseEinstein(Eqm_Temp, omega2) -(DoS.gox_orig[freqID]*BoseEinstein(Eqm_Temp, w)) #(GoX_omega1+GoX_omega2-DoS.gox_orig[freqID])*(BoseEinstein(Eqm_Temp,omega2)+DoS.gox[w1ID]-BoseEinstein(Eqm_Temp,w))

                   # num_doesexist = (GoX_omega1*GoX_omega2)*(BoseEinstein(ShTemp,omega2)*BoseEinstein(Eqm_Temp,omega1)-BoseEinstein(Eqm_Temp,w)) #-DoS.gox[freqID] #(BoseEinstein(Eqm_Temp,omega1)*BoseEinstein(Eqm_Temp,omega2)-BoseEinstein(Eqm_Temp,DoS.frequency[freqID])) #*BoseEinstein(Eqm_Temp,w))
##THESE TWO LINES ARE EXPERIMENTAL
                   try:
                       REVERSE
                   except NameError:
                       REVERSE=1
                   if REVERSE==0:
                      if num_doesexist < 0:
                          num_doesexist =  0
                       
                else:
                    num_doesexist=0
          #This next section is if you want to pump above 3Wmax from scattering between a phonon and a 2-3Wmax state..          
#        elif w >= 2*WMax and omega2 <= WMax and omega1 >= 2*WMax and omega1<w: # PHONONConserved==1: ### omega2 <= omega1 OR omega2 < omega1
#               # print('made it! w=' , w , ' W1= ', omega1)
#                w2IDint = np.where(np.isclose(DoS.frequency,omega2,atol=toler))
#                
#                if isinstance(w2IDint, int) == False:
#                    w2IDint=w2IDint[0][0]  #This gives us the element position of omega2 in the frequency array
#                    
#                w1IDint=np.where(DoS.frequency==omega1)
#                w1ID=w1IDint[0][0] 
#                
#                w2ID=w2IDint
#                GoX_omega2 = DoS.gox_orig[w2ID]
#                GoX_omega1 = DoS.gox_orig[w1ID]
#               
#                if GoX_omega2 > intcut and GoX_omega1 > intcut:
##                else:
#               # if GoX_omega1 > intcut and GoX_omega2 > intcut and DoS.gox[freqID] > intcut:
#                    num_doesexist = (GoX_omega1+GoX_omega2-DoS.gox_orig[freqID])*(BoseEinstein(Eqm_Temp,omega2)+BoseEinstein(Eqm_Temp,omega1)-BoseEinstein(Eqm_Temp,w))
#                   # num_doesexist = (GoX_omega1*GoX_omega2)*(BoseEinstein(ShTemp,omega2)*BoseEinstein(Eqm_Temp,omega1)-BoseEinstein(Eqm_Temp,w)) #-DoS.gox[freqID] #(BoseEinstein(Eqm_Temp,omega1)*BoseEinstein(Eqm_Temp,omega2)-BoseEinstein(Eqm_Temp,DoS.frequency[freqID])) #*BoseEinstein(Eqm_Temp,w))
#                else:
#                    num_doesexist=0
      
                if omega1==omega2 and CHECK==0:
                    print('-------------------')
                    print('Omega1: ', omega1)
                    print('Omega2: ', omega2)
                    print('numdoesexist:  ', num_doesexist)
                    print('GoX_w1:  ', GoX_omega1)
                    print('GoX_w2:  ', GoX_omega2)
                    print('none:  ', n_one)
                    print('ntwo:  ', n_two)
                    print('-------------------')
                    #print('updating values 2pdos', DoS.tpd_temp)
        else: 
            num_doesexist = 0 # print('Does it change -before   ' , DoS.tpd_temp)
        
        DoS.tpd_temp +=num_doesexist  #line 324 in matlab

    DoS.twoPDoS[freqID]= DoS.tpd_temp

try:
    otone_num
except NameError:
    otone_num=1
    print('Set Otone_Num to 1 - seems you did not run the two-step model')

Proj_DoS(DoS,2,otone_num,PBVal)

#### INTEGRATE THE Projected 2pDOS
MinIntegration=float(input('Minimum integration N*WMax: '))
#MaxIntegrationID = np.where(DoS.twoPDoS[np.int(2*WMax):]<intcut)[0][0]
#MaxIntegration = DoS.frequency[np.int(2*WMax)+MaxIntegrationID]    
MaxIntegration=float(input('Maximum integration N*WMax: '))
IntIDBottom = np.where(np.isclose(DoS.frequency,MinIntegration*WMax,atol=toler))
IntIDBottom = IntIDBottom[0][0]

#IntIDTop = np.where(np.isclose(DoS.frequency,MaxIntegration*1,atol=toler))
IntIDTop = np.where(np.isclose(DoS.frequency,MaxIntegration*WMax,atol=toler))
if isinstance(IntIDTop, int) == False:
    IntIDTop = IntIDTop[0][0]
#if MaxIntegration*WMax > np.max(DoS.frequency):
#    IntIDTop = np.size(DoS.frequency)-1
    

print('Bottom integration frequency:  ', DoS.frequency[IntIDBottom])
print('Top integration frequency:  ', DoS.frequency[IntIDTop])  

try:
    print_norm
except NameError:
    print_norm=1

if print_norm ==0:
    int_projtwopdos=integrate.trapz(DoS.proj_twoPDoS[IntIDBottom:IntIDTop])   ##Check exactly which index these come in in Python conventions
    print('Non-Normalised Integration of 2pDOS is:  ', int_projtwopdos)
    WNorm_int_projtwopdos = int_projtwopdos / (DoS.frequency[IntIDTop]-DoS.frequency[IntIDBottom])
    print('DoS-Wrange-Normalised Integration of 2pDOS is:  ', WNorm_int_projtwopdos)
    Norm_int_projtwopdos=int_projtwopdos / integrate.trapz(DoS.gox_orig)
    print('DoS-Normalised Integration of 2pDOS is:  ', Norm_int_projtwopdos)
    NormNorm_int_projtwopdos = Norm_int_projtwopdos / (DoS.frequency[IntIDTop]-DoS.frequency[IntIDBottom])  
    print('DoS-Wrange-Normalised-3N-Normalised of 2pDOS is:  ', NormNorm_int_projtwopdos)



#MatLabInt = np.trapz(DoS.proj_twoPDoS[IntIDBottom:IntIDTop])/(np.trapz(DoS.gox)*np.trapz(DoS.gox[IntIDBottom:IntIDTop]))
MatLabInt=integrate.simps(DoS.proj_twoPDoS[IntIDBottom:IntIDTop]) #10**(6) #/(1*np.trapz(DoS.gox_orig[IntIDBottom:IntIDTop]))*10**(-4)
print('Integration over portion of states in up-pumped region:  ', MatLabInt)


NBOND=np.int(input('Number of explosophoric bonds per molecule:  '))
MatLabInt2=integrate.simps(DoS.proj_twoPDoS[IntIDBottom:IntIDTop])/NBOND
print('Integrated normalised to number of explosophoric bonds:  ', MatLabInt2)

TOTAL_ENERGY=NtoEnergy(DoS)
IntTotal_Energy = integrate.simps(DoS.VibEnergy[IntIDBottom:IntIDTop])
print('total energy in region : ', IntTotal_Energy)

###PLOTTING###
## note will write on top of each other, so comment out what you don't want
plt.figure('Two Phonon DOS')
plt.plot(DoS.frequency, DoS.twoPDoS,linestyle='-',label='TwoPDoS Wmax=%.1f'%WMax)
plt.legend()
plt.figure('Projected DOS')
plt.plot(DoS.frequency, DoS.proj_twoPDoS,linestyle='--',label='Proj2pDOS Wmax=%.1f'%WMax)
plt.plot(DoS.otone_freq,DoS.proj_otone,linestyle='-',label='ProjOtone Wmax=%.1f'%WMax)
plt.legend()
plt.figure('Projected Energy')
plt.plot(DoS.frequency,DoS.VibEnergy,linestyle='-',label='Energy Wmax=%.1f'%WMax)


###WRITE OUTPUT TO CSV###
if 'y' in input('save files?:  '):
    data=[DoS.frequency, DoS.gox_orig, DoS.gox, DoS.twoPDoS, DoS.proj_twoPDoS, DoS.otone_freq, DoS.proj_otonein,DoS.ShockAll ]

    main_window = Tk()

    a = asksaveasfilename(filetypes=(("PNG Image", "*.png"),("All Files", "*.*")), 
                          defaultextension='.png', title="Save Up-pumping figure")
    b = asksaveasfilename(filetypes=(("CSV", "*.csv"),("All Files", "*.*")), 
                          defaultextension='.csv', title="Save Up-pumping data")
    if a:
        plt.savefig(a)
        if b:
            out = csv.writer(open(b,"w", newline=''), delimiter=',')
            out.writerows(zip(*data))
            out=csv.writer
#ff=asksaveasfile(initialdir=wd, filetypes=(('png files', '*.png'), ('jpg files', '*.jpg')),title="Window-2")
#if ff:
#    plt.savefig(ff)
    main_window.destroy()
