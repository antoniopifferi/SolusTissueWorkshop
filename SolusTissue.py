# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 09:18:29 2021

@author: anton
"""
from numpy import *
#from scipy import *
import matplotlib
#matplotlib.rcParams['text.usetex'] = True
from matplotlib.pyplot import *
from pandas import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
close("all")

# PARAMETERS
PATH_MAIN='C:\\OneDrivePolimi\\OneDrive - Politecnico di Milano\\Beta\\Analysis\\SolusTissue\\'
PATH_DATA='Data\\'
PATH_RESULTS='Results\\'
PATH_SETTINGS='Settings\\'

FILE_DATA='Data1.txt'
FILE_LABBOOK='Labbook1.txt'
FILE_SCENARIO='Scenario1.txt'
FILE_VARIABLE='Variable1.txt'

# FILE_DATA='_Meat_FoM_spectral_reconcurves_fromSpectralFit2_.txt'
# FILE_LABBOOK='LabbookMeat3.txt'
# FILE_SCENARIO='ScenarioMeat3.txt'
# FILE_VARIABLE='VariableMeat3.txt'

FIGWIDTH=15
SAVE_FIG=False
SUP_TITLE=True
ASPECT_RATIO=True

# LOAD
Data=read_csv(PATH_MAIN+PATH_DATA+FILE_DATA,sep='\t')
Labbook=read_csv(PATH_MAIN+PATH_SETTINGS+FILE_LABBOOK,sep='\t')
Scenario=read_csv(PATH_MAIN+PATH_SETTINGS+FILE_SCENARIO,sep='\t')  
Variable=read_csv(PATH_MAIN+PATH_SETTINGS+FILE_VARIABLE,sep='\t')  

# COMPLETE DATA
Variable.Label.fillna(Variable.NewVar,inplace=True)
Variable.Unit.fillna("",inplace=True)
dcVariable=dict(zip(Variable.OldVar, Variable.NewVar))
dcUnit=dict(zip(Variable.NewVar, Variable.Unit))
dcLabel=dict(zip(Variable.NewVar, Variable.Label))
Data.rename(columns=dcVariable,inplace=True)
Data = Data.merge(Labbook, on='Device')

# CALC VARIABLES
#Data['bkgMuaTrue']=1/(1+Data['contrMuaTrue'])*Data['incMuaTrue'] # TRUE BKG MUA
for var,fact in zip(Variable[Variable.Factor>0].NewVar,Variable[Variable.Factor>0].Factor): Data[var]=Data[var]*fact

# PLOT
for i,s in Scenario.iterrows(): # iterate over the whole Scenario
    
    # extract arrays
    if notnull(s.Truth): Data[s.Truth+'1']=Data[s.Truth]
    DataE=Data
    if notnull(s.Extract1): DataE=DataE[DataE[s.Extract1]==s.eVal1]
    aCol=DataE[s.Col].unique()
    aRow=DataE[s.Row].unique()    
    Name=s.Var+"_"+s.View


    #do plot
    nRow=len(aRow)
    nCol=len(aCol)
    aRatio = (nRow+0.5)/(nCol+1) if ASPECT_RATIO else 9/16
    figwidth = FIGWIDTH*0.6 if nCol==3 else FIGWIDTH
    figData,axs=subplots(nRow,nCol,num='Fig'+str(Name),figsize=(figwidth,aRatio*figwidth),squeeze=False)
    if SUP_TITLE: suptitle(FILE_SCENARIO+'  #  '+FILE_DATA+'  #  '+str(Name))
    for iCol,oCol in enumerate(aCol):
        for iRow,oRow in enumerate(aRow):
            axi=axs[iRow,iCol]
            sca(axi)

            subData=DataE[(DataE[s.Col]==oCol)&(DataE[s.Row]==oRow)]
            table=subData.pivot_table(values=s.Y,index=s.X,columns=s.Line,aggfunc='mean')
            table.style.format({'PertMua':'{0:,.0f} nm','horsepower':'{0:,.0f}hp'})
            table.plot(ax=axi,marker='D',legend=False,xlabel=False)
            if((iCol==0)and(iRow==0)): legend()
            if notnull(s.Truth): truth=subData.pivot_table(values=s.Truth+'1',index=s.X,columns=s.Line,aggfunc='mean')
            if notnull(s.Truth): truth.plot(ax=axi,marker='',color='black',legend=False)
            if notnull(s.Ymin): gca().set_ylim([s.Ymin,s.Ymax]) # check if there is any value, including 0, otherwise leave autoscale
            grid(True)

            # print labels
            xLab = (dcLabel[s.X] if s.X in dcLabel else s.X) + (" ("+dcUnit[s.X]+")" if s.X in dcUnit else "")
            yLab = (dcLabel[s.Y] if s.Y in dcLabel else s.Y) + (" ("+dcUnit[s.Y]+")" if s.Y in dcUnit else "")
            rLab = (dcLabel[s.Row] if s.Row in dcLabel else s.Row)+"="+str(oRow) + (" "+dcUnit[s.Row] if s.Row in dcUnit else "")
            cLab = (dcLabel[s.Col] if s.Col in dcLabel else s.Col)+"="+str(oCol) + (" "+dcUnit[s.Col] if s.Col in dcUnit else "")    
            if iCol==0: gca().set_ylabel(yLab)
            if iRow==(nRow-1): gca().set_xlabel(xLab)           
            if iRow==0: gca().set_title(cLab)
            if iCol==(nCol-1): gca().twinx().set_ylabel(rLab)
            
            
    figData.tight_layout()
    show()
    if SAVE_FIG: figData.savefig(PATH_D48+PATH_RESULTS+'Fig_'+str(Name)+'.jpg',format='jpg')
