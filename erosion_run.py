# -*- coding: utf-8 -*-
import numpy as np
from netCDF4 import Dataset
import sys
sys.path.insert(0,'C:\VIKA\Postdoc_LSCE\Indices\E_code')
import erosion

###########################################
#This script calculates global yearly average soil erosion rates
#Transient simulation = changing land use and climate (CC+LUC), changing climate only (CC),
#or changing land use only (LUC)
#Equilibrium simulation = average soil erosion rates for a specific period
#When using this model code please cite: "Naipal ,V., Ciais, P., Wang, Y., Lauerwald, R.,
#Guenet, B., Van Oost, K. (2018). Global soil organic carbon removal by water erosion under
#climate change and land use change during 1850–2005 AD. Biogeosciences, 15, 4459–4480,
#https://doi.org/10.5194/bg-15-4459-2018"
############################################

print ('start script')

#settings
workdir="C:\VIKA\Postdoc_LSCE\Indices" #Enter here your work directory
years=range(1975,2006) #transient period
simulation="CC+LUC" #select simulation type from: equilibrium,CC+LUC,CC,LUC
clim_class=("A","B","C","D")
clim_sel=1 #select climate classification: Zone A, B, C or D. Choose from: 1, 2, 3 or 4
lc=1 #select vegetation type crop, 0=tree, 2=grass, 3=bare
veget="crop" #select vegetation type tree,crop,grass or bare

rows=1680;kols=4322 #global 5 arcmin resolution grid

#general input data
data = Dataset('%s\input\climate_classification_aggr_5m.nc' % (workdir) ,'r') #climate zone data
clim= data.variables['climate'][:] #1= Zone A (tropical), 2= Zone B (arid/dry), 3= Zone C (temperate), 4= Zone D (Cold)
clim[clim<0.]=np.nan
clim = np.ravel(clim)
data.close()  

data = Dataset('%s\input\S_scaled_5m.nc' % (workdir),'r') #S factor
S= data.variables['s'][:]
S[S<0.] = 0.
S[np.isnan(S)==True] = 0.
S=np.ravel(S)
S[clim!=clim_sel]=0.
data.close()

data = Dataset('%s\input\K_volcanic_nogravel_5m.nc' % (workdir),'r')
K = data.variables['k'][:] # K factor without gravel correction but with volcanic soil correction
K[K<0.] = 0.
K[np.isnan(K)==True] = 0.
K=np.ravel(K)
K[clim!=clim_sel]=0.
data.close()

data = Dataset('%s\input\gravel_topsoil_5m.nc' % (workdir),'r')
gravel = data.variables['GRAV'][:] # gravel content in %
gravel=gravel.astype(float)
gravel[gravel<0.] = 0.
gravel[np.isnan(gravel)==True] = 0.
gravel = np.ravel(gravel)
gravel[clim!=clim_sel]=0.
data.close()

print ('read input data part 1: done')

if simulation =='equilibrium':
  ####################################
  # Equilibrium 
  ###################################
  data0 = Dataset('%s\input\landcover\PFTmap_LUHv2_BM3_HoughtonCountryForestarea_crop_1975-1985_mean_5m.nc' % (workdir),'r') #crop fraction
  crop = data0.variables['maxvegetfrac'][0]
  crop[crop>1.] = np.nan;crop[crop<0.] = np.nan
  crop=np.nansum(crop,axis=0)
  crop[np.isnan(crop)==True] = 0.
  crop=np.ravel(crop)
  crop[clim!=clim_sel]=0.
  data0.close()

  data2 = Dataset('%s\input\R_factor\ISIMIP2b_R_1975-1985_mean_5m.nc' % (workdir),'r') # R factor
  R = data2.variables['r'][:] 
  R[R<0]=0.;R[R>10e35]=0.;R[np.isnan(R)==True] = 0.
  R=np.ravel(R)
  R[clim!=clim_sel]=0.
  data2.close()    

  nc = Dataset('%s\input\C_factor\C_crop_1975-1985_mean_5m.nc' % (workdir),'r') # C factor   
  C = nc.variables['c_lc'][:]
  C[C<0]=0.;C[np.isnan(C)==True] = 0.
  C=np.ravel(C)   
  C[clim!=clim_sel]=0.
  
  print ('read input data part 2: done')
  
  #calculate erosion for this timeperiod
  E=erosion.erosion.main(S,K,R,C,veget,crop,gravel,rows,kols)
  E = E.reshape(rows,kols)      
  #output
  output = Dataset('%s\output\erosion\E_%s_mean_1975-1985_zone%s.nc' % (workdir,veget,clim_class[clim_sel-1]),'w')
  output.createDimension('latitude',rows)
  output.createDimension('longitude',kols)
  var=output.createVariable('latitude','f',('latitude',))
  var[:]=[83.9917-n*0.08329773 for n in range(rows)]
  var=output.createVariable('longitude','f',('longitude',))
  var[:]=[-180.+n*0.0833 for n in range(kols)]
  output.createVariable('E','d',('latitude','longitude',))
  output.variables['E'][:] = E
  output.close()  

  print ('erosion output: done')
    
  nc.close()    

if simulation =='CC+LUC':
  ####################################
  # CC + LUC
  ###################################
  for t in years:
    data = Dataset('%s/input/R_factor/ISIMIP2b_R_%04i_5m.nc' % (workdir,t),'r') # R factor
    R = data.variables['r'][:] 
    R[R<0]=0.;R[R>10e35]=0.
    R[np.isnan(R)==True] = 0.
    R=np.ravel(R)
    data.close()    

    data = Dataset('%s/input/landcover/PFTmap_LUHv2_BM3_HoughtonCountryForestarea_crop_%04i_5m.nc' % (workdir,t),'r')  #crop fraction
    crop = data.variables['maxvegetfrac'][0]
    crop[crop>1.] = np.nan;crop[crop<0.] = np.nan
    crop=np.nansum(crop,axis=0)
    crop[np.isnan(crop)==True] = 0.
    crop=np.ravel(crop)
      
    data = Dataset('%s/input/C_factor/C_crop_%04i_5m.nc' % (workdir,t),'r')
    C = data.variables['c_lc'][:] # C factor for crop
    C[C<0]=0.
    C[np.isnan(C)==True] = 0.
    C=np.ravel(C)   
    data.close()
    
    #calculate erosion for this timeperiod
    E=erosion.erosion.main(S,K,R,C,veget,crop,gravel,rows,kols)
    E = E.reshape(rows,kols)
    
    #output 
    output = Dataset('%s/output/erosion/E_%s_%04i_zone%s.nc' % (workdir,veget,t,clim_class[clim_sel-1]),'w')
    output.createDimension('latitude',rows)
    output.createDimension('longitude',kols)
    var=output.createVariable('latitude','f',('latitude',))
    var[:]=[83.9917-n*0.08329773 for n in range(rows)]
    var=output.createVariable('longitude','f',('longitude',))
    var[:]=[-180.+n*0.0833 for n in range(kols)]
    output.createVariable('E','d',('latitude','longitude',))
    output.variables['E'][:] = E
    output.close() 

if simulation =='CC':
  ####################################
  # CC only (LUC=set to 1975-1985)
  ###################################
  data0 = Dataset('%s/input/landcover/PFTmap_LUHv2_BM3_HoughtonCountryForestarea_crop_1975-1985_mean_5m.nc' % (workdir),'r')  #crop fraction
  crop = data0.variables['maxvegetfrac'][0]
  crop[crop>1.] = np.nan;crop[crop<0.] = np.nan
  crop=np.nansum(crop,axis=0)
  crop[np.isnan(crop)==True] = 0.
  crop=np.ravel(crop)
  data0.close()

  for t in years:
    data2 = Dataset('%s/input/R_factor/ISIMIP2b_R_%04i_5m.nc' % (workdir,t),'r') # R factor
    R = data2.variables['r'][:] 
    R[R<0]=0.
    R[R>10e35]=0.
    R[np.isnan(R)==True] = 0.
    R=np.ravel(R)
    data2.close()    
      
    data1 = Dataset('%s/input/C_factor/CC/C_crop_%04i_5m.nc' % (workdir,t),'r') # C factor
    C = data1.variables['c_lc'][:]
    C[C<0]=0.;C[np.isnan(C)==True] = 0.
    C=np.ravel(C)   
      
    #calculate erosion for this timeperiod
    E=erosion.erosion.main(S,K,R,C,veget,crop,gravel,rows,kols)
    E = E.reshape(rows,kols)      
    #output
    output = Dataset('%s/output/erosion/CC/E_%s_%04i_zone%s.nc' % (workdir,veget,t,clim_class[clim_sel-1]),'w')
    output.createDimension('latitude',rows)
    output.createDimension('longitude',kols)
    var=output.createVariable('latitude','f',('latitude',))
    var[:]=[83.9917-n*0.08329773 for n in range(rows)]
    var=output.createVariable('longitude','f',('longitude',))
    var[:]=[-180.+n*0.0833 for n in range(kols)]
    output.createVariable('E','d',('latitude','longitude',))
    output.variables['E'][:] = E
    output.close()   

if simulation =='LUC':
  ####################################
  # LUC only (CC=set to 1975-1985)
  ###################################
  data0 = Dataset('%s/input/R_factor/ISIMIP2b_R_1975-1985_mean_5m.nc' % (workdir),'r') # R factor
  R = data0.variables['r'][:] 
  R[R<0]=0.;R[R>10e35]=0.
  R[np.isnan(R)==True] = 0.
  R=np.ravel(R)
  data0.close() 

  for t in years: 
    data1 = Dataset('%s/input/landcover/PFTmap_LUHv2_BM3_HoughtonCountryForestarea_crop_%04i_5m.nc' % (workdir,t),'r') #crop fraction
    crop = data1.variables['maxvegetfrac'][0]
    crop[crop>1.] = np.nan;crop[crop<0.] = np.nan
    crop=np.nansum(crop,axis=0)
    crop[np.isnan(crop)==True] = 0.
    crop=np.ravel(crop)
    data1.close()
    
    data1 = Dataset('%s/input/C_factor/LUC/C_crop_%04i_5m.nc' % (workdir,t),'r') # C factor
    C = data1.variables['c_lc'][:]
    C[C<0]=0.;C[np.isnan(C)==True] = 0.
    C=np.ravel(C)   
      
    #calculate erosion for this timeperiod
    E=erosion.erosion.main(S,K,R,C,veget,crop,gravel,rows,kols)
    E = E.reshape(rows,kols)      
    #output
    output = Dataset('%s/output/erosion/LUC/E_%s_%04i_zone%s.nc' % (workdir,veget,t,clim_class[clim_sel-1]),'w')
    output.createDimension('latitude',rows)
    output.createDimension('longitude',kols)
    var=output.createVariable('latitude','f',('latitude',))
    var[:]=[83.9917-n*0.08329773 for n in range(rows)]
    var=output.createVariable('longitude','f',('longitude',))
    var[:]=[-180.+n*0.0833 for n in range(kols)]
    output.createVariable('E','d',('latitude','longitude',))
    output.variables['E'][:] = E
    output.close()        

