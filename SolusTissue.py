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
PATH_MAIN='C:\\OneDrivePolimi\\OneDrive - Politecnico di Milano\\Beta\\Analysis\\SolusTissueWorkshop\\'
PATH_DATA='Data\\'
PATH_RESULTS='Results\\'
PATH_SETTINGS='Settings\\'

FILE_DATA='Data1.txt'
FILE_LABBOOK='Labbook1.txt'
FILE_SCENARIO='Scenario3.txt'
FILE_VARIABLE='Variable1.txt'
FILE_COMPONENTS='Components1.txt'
FILE_CONCENTRATIONS='ConcentrationsSOLUS1.txt'
DEVICE_CONCENTRATIONS='SOLUS'

# FILE_DATA='_Meat_FoM_spectral_reconcurves_fromSpectralFit2_.txt'
# FILE_LABBOOK='LabbookMeat3.txt'
# FILE_SCENARIO='ScenarioMeat3.txt'
# FILE_VARIABLE='VariableMeat3.txt'

#AcceptVolunteer=[1,4,6,7,8,9,10]
#AcceptTissue=['Forehead','Deltoid','Calcaneus','Abdomen']
AcceptVolunteer=[7,8,10]
AcceptTissue=['Forehead','Deltoid','Abdomen','Calcaneus']
AcceptRho=[2.6]
AcceptLambda=arange(630,1100,10)

FIGWIDTH=15
FONT_SIZE=13
SAVE_FIG=True
SUP_TITLE=False
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

# FILTER DATA
#Data=Data[Data.Volunteer.isin(AcceptVolunteer)]
#Data=Data[Data.Tissue.isin(AcceptTissue)]
Data=Data[Data.Rho.isin(AcceptRho)]
Data=Data[Data.Lambda.isin(AcceptLambda)]


# CALC VARIABLES
#Data['bkgMuaTrue']=1/(1+Data['contrMuaTrue'])*Data['incMuaTrue'] # TRUE BKG MUA
for var,fact in zip(Variable[Variable.Factor>0].NewVar,Variable[Variable.Factor>0].Factor): Data[var]=Data[var]*fact

# PLOT
matplotlib.rcParams.update({'font.size': FONT_SIZE})
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
    # figData,axs=subplots(nRow,nCol,num='Fig'+str(Name),figsize=(figwidth,aRatio*figwidth),squeeze=False)
    figData,axs=subplots(nRow,nCol,num='Fig'+str(Name),sharex=True,sharey=True,figsize=(figwidth,aRatio*figwidth),squeeze=False)
    subplots_adjust(hspace=0,wspace=0,left=0.09,bottom=0.09)
#    figData.supxlabel(fontsize=12,fontweight="bold")
    # figData.supxlabel(x_axis[s]+" ("+units[x_axis[s]]+")",fontsize=12,fontweight="bold")
# fig.supylabel(y_axis[s]+" ("+units[y_axis[s]]+")",fontsize=12,fontweight="bold")
# fig.suptitle("%s"%title[s]+"  %s"%page[s]+"=%s"%p,fontsize=12,fontweight = 'bold')            
    if SUP_TITLE: suptitle(FILE_SCENARIO+'  #  '+FILE_DATA+'  #  '+str(Name))


# fig,axs = plt.subplots(nrows,ncols,sharex=True,sharey=False,figsize=(width,hieght),squeeze=False)
# plt.subplots_adjust(hspace=0,wspace=0,left=0.09,bottom=0.09)
# fig.supxlabel(x_axis[s]+" ("+units[x_axis[s]]+")",fontsize=12,fontweight="bold")
# fig.supylabel(y_axis[s]+" ("+units[y_axis[s]]+")",fontsize=12,fontweight="bold")
# fig.suptitle("%s"%title[s]+"  %s"%page[s]+"=%s"%p,fontsize=12,fontweight = 'bold')            




    for iCol,oCol in enumerate(aCol):
        for iRow,oRow in enumerate(aRow):
            axi=axs[iRow,iCol]
            sca(axi)

            subData=DataE[(DataE[s.Col]==oCol)&(DataE[s.Row]==oRow)]
            table=subData.pivot_table(values=s.Y,index=s.X,columns=s.Line,aggfunc='mean')
            #table.style.format({'PertMua':'{0:,.0f} nm','horsepower':'{0:,.0f}hp'})
            #table.plot(ax=axi,marker='D',legend=False,xlabel=False)
            table.plot(ax=axi,marker='D',legend=False)
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
            # if iCol==0: gca().set_ylabel(yLab)
            # if iRow==(nRow-1): gca().set_xlabel(xLab)           
            if iRow==0: gca().set_title(cLab)
            if iCol==(nCol-1): gca().twinx().set_ylabel(rLab)
            #gca().set_ylabel(yLab)
            #gca().set_xlabel(xLab)
 
    # SAVE FIGURE        
    figData.tight_layout()
    show()
    if SAVE_FIG: figData.savefig(PATH_MAIN+PATH_RESULTS+'Fig_'+str(Name)+'.jpg',format='jpg')


#%% CALC COMPONENTS
Components=read_table(PATH_MAIN+PATH_DATA+FILE_COMPONENTS)
figure(num='FigComp')
Components.plot(x='Lambda')
yscale('log')
ylim([0,0.5]), title('Components'), xlabel('wavelength (nm)'), ylabel('specific absorption (cm-1)')
show()

# gbData=Data.groupby(['Device','Tissue','Volunteer','Lambda'])['Mua']
# meanVal=gbData.mean()
table=Data[Data.Device==DEVICE_CONCENTRATIONS].pivot_table(values='Mua',index='Lambda',columns=['Device','Tissue','Volunteer'],aggfunc='mean')
comp=Components[Components['Lambda'].isin(Data[Data.Device==DEVICE_CONCENTRATIONS].Lambda.unique())].values
comp=delete(comp,0,1)
aComp=linalg.lstsq(comp,table,rcond=None)[0] #[0] to extract m-coeff in lstsq
# y=log(table.loc[LAMBDA1:LAMBDA2,'Mus'])
#y=log(table.loc[LAMBDA1:LAMBDA2,'Mus'])
# x=log(y.index/LAMBDA0)
# model=polyfit(x,y,1)
# b=-model[0]
# a=exp(model[1])
# #A = vstack([x, np.ones(len(x))]).T
dfComp=DataFrame(data=aComp.transpose(),columns=Components.columns[1:],index=table.columns)
dfComp['tHb']=dfComp['HHb']+dfComp['O2Hb']
dfComp['SO2']=dfComp['O2Hb']/dfComp['tHb']
dfComp['Tot']=dfComp['Lipid']+dfComp['H2O']
# dfComp['Tot']=dfComp['Lipid']+dfComp['H2O']+dfComp['Coll']
dfComp['FitComp']='LambdaFit'
# dfComp['a']=a
# dfComp['b']=b
# dfComp.plot()
dfComp.to_csv(path_or_buf=PATH_MAIN+PATH_RESULTS+FILE_CONCENTRATIONS,sep='\t')
# figure(num='FigConc')
# plot(x,y)
# #filtData=merge(filtData,dfComp,on=['Subject','Meas','Rho'])

# show()




# # SINGLE FIGURE
# #figTissue=figure(num='FigTissue',figsize=(figwidth,aRatio*figwidth))
# figTissue=figure(num='FigTissue')
# #filtData=Data[(Data.Volunteer==Volunteer0) && (Data.Tissue==selectTissue)]
# table=Data[(Data.Device=='DOS') & (Data.Volunteer==9)].pivot_table(values='Mua',index='Lambda',columns='Tissue',aggfunc='mean')
# ax0=table.plot()
# table=Data[(Data.Device=='SOLUS') & (Data.Volunteer==9)].pivot_table(values='Mua',index='Lambda',columns='Tissue',aggfunc='mean')
# gca().set_prop_cycle(None)
# table.plot(ax=ax0,marker='D',linestyle='None',legend=False)
# ylim([0,0.5])
# if SAVE_FIG: savefig(PATH_MAIN+PATH_RESULTS+'Fig_'+str(Name)+'.jpg',format='jpg')

# figVol=figure(num='FigVol')
# #filtData=Data[(Data.Volunteer==Volunteer0) && (Data.Tissue==selectTissue)]
# table=Data[(Data.Device=='DOS') & (Data.Tissue=='Calcaneus')].pivot_table(values='Mua',index='Lambda',columns='Volunteer',aggfunc='mean')
# ax0=table.plot()
# table=Data[(Data.Device=='SOLUS') & (Data.Tissue=='Calcaneus')].pivot_table(values='Mua',index='Lambda',columns='Volunteer',aggfunc='mean')
# gca().set_prop_cycle(None)
# table.plot(ax=ax0,marker='D',linestyle='None',legend=False)
# ylim([0,0.5])
# if SAVE_FIG: savefig(PATH_MAIN+PATH_RESULTS+'Fig_'+str(Name)+'.jpg',format='jpg')


