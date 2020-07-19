# -*- coding: utf-8 -*-
import numpy as np
from netCDF4 import Dataset

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
workdir= '/home/wieka20/Dokumente/Postdoc_LMU/Teaching/Practical/' #Enter here your work directory
years=range(1975,2006) #transient period
reg='EU_west' #select your region, 'EU_east' for Eastern Europe, 'EU_west' for Western Europe, 'BRA' for Brazil or 'AUS' for Australia
simulation="CC+LUC" #select simulation type from: equilibrium,CC+LUC,CC,LUC
lc=1 #select vegetation type crop, 0=tree, 2=grass, 3=bare
veget="crop" #select vegetation type tree,crop,grass or bare

#general input data
data = Dataset('%s/input_%s/rest/S_scaled_5m_%s.nc' % (workdir,reg,reg),'r') #S factor
S= data.variables['s'][:]
lat0=data.variables['lat'][0]
lon0=data.variables['lon'][0]
rows=len(S[:,0]);kols=len(S[0]) #global 5 arcmin resolution grid
S[S<0.] = 0.
S[np.isnan(S)==True] = 0.
S=np.ravel(S)
data.close()

data = Dataset('%s/input_%s/rest/K_volcanic_nogravel_5m_%s.nc' % (workdir,reg,reg),'r')
K = data.variables['k'][:] # K factor without gravel correction but with volcanic soil correction
K[K<0.] = 0.
K[np.isnan(K)==True] = 0.
K=np.ravel(K)
data.close()

data = Dataset('%s/input_%s/rest/gravel_topsoil_5m_%s.nc' % (workdir,reg,reg),'r')
gravel = data.variables['GRAV'][:] # gravel content in %
gravel=gravel.astype(float)
gravel[gravel<0.] = 0.
gravel[np.isnan(gravel)==True] = 0.
gravel = np.ravel(gravel)
data.close()

print ('read input data part 1: done')

if simulation =='equilibrium':
  ####################################
  # Equilibrium 
  ###################################

  data2 = Dataset('%s\input_%s\R_factor\ISIMIP2b_R_1975-1985_mean_5m_%s.nc' % (workdir,reg,reg),'r') # R factor
  R = data2.variables['r'][:] 
  R[R<0]=0.;R[R>10e35]=0.;R[np.isnan(R)==True] = 0.
  R=np.ravel(R)
  data2.close()    

  nc = Dataset('%s\input_%s\C_factor\C_crop_1975-1985_mean_5m_%s.nc' % (workdir,reg,reg),'r') # C factor   
  C = nc.variables['c_lc'][:]
  C[C<0]=0.;C[np.isnan(C)==True] = 0.
  C=np.ravel(C)   
  
  print ('read input data part 2: done')
  
  #calculate erosion for this timeperiod
  E=S*R*K*C #potential erosion
  E[gravel>12.]=E-(0.8*E) #erosion reduction by 80% if gravel on cropland >12%
  E[E>100.]=100. #max erosion rate is 100t/ha/y
  E = E.reshape(rows,kols)      
  #output
  output = Dataset('%s\output\erosion\E_%s_mean_1975-1985_%s.nc' % (workdir,veget,reg),'w')
  output.createDimension('latitude',rows)
  output.createDimension('longitude',kols)
  var=output.createVariable('latitude','f',('latitude',))
  var[:]=[lat0-n*0.0833 for n in range(rows)]
  var=output.createVariable('longitude','f',('longitude',))
  var[:]=[lon0+n*0.0833 for n in range(kols)]
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
    data = Dataset('%s/input_%s/R_factor/ISIMIP2b_R_%04i_5m_%s.nc' % (workdir,reg,t,reg),'r') # R factor
    R = data.variables['r'][:] 
    R[R<0]=0.;R[R>10e35]=0.
    R[np.isnan(R)==True] = 0.
    R=np.ravel(R)
    data.close()    
      
    data = Dataset('%s/input_%s/C_factor/C_crop_%04i_5m_%s.nc' % (workdir,reg,t,reg),'r')
    C = data.variables['c_lc'][:] # C factor for crop
    C[C<0]=0.
    C[np.isnan(C)==True] = 0.
    C=np.ravel(C)   
    data.close()
    
    #calculate erosion for this timeperiod
    E=S*R*K*C #potential erosion
    E[gravel>12.]=E[gravel>12.]-(0.8*E[gravel>12.]) #erosion reduction by 80% if gravel on cropland >12%
    E[E>100.]=100. #max erosion rate is 100t/ha/y
    E = E.reshape(rows,kols) 
    
    #output 
    output = Dataset('%s/output/erosion/E_%s_%04i_%s.nc' % (workdir,veget,t,reg),'w')
    output.createDimension('latitude',rows)
    output.createDimension('longitude',kols)
    var=output.createVariable('latitude','f',('latitude',))
    var[:]=[lat0-n*0.0833 for n in range(rows)]
    var=output.createVariable('longitude','f',('longitude',))
    var[:]=[lon0+n*0.0833 for n in range(kols)]
    output.createVariable('E','d',('latitude','longitude',))
    output.variables['E'][:] = E
    output.close() 

if simulation =='CC':
  ####################################
  # CC only (LUC=set to 1975-1985)
  ###################################
  for t in years:
    data2 = Dataset('%s/input_%s/R_factor/ISIMIP2b_R_%04i_5m_%s.nc' % (workdir,reg,t,reg),'r') # R factor
    R = data2.variables['r'][:] 
    R[R<0]=0.
    R[R>10e35]=0.
    R[np.isnan(R)==True] = 0.
    R=np.ravel(R)
    data2.close()    
      
    data1 = Dataset('%s/input_%s/C_factor/CC/C_crop_%04i_5m_%s.nc' % (workdir,reg,t,reg),'r') # C factor
    C = data1.variables['c_lc'][:]
    C[C<0]=0.;C[np.isnan(C)==True] = 0.
    C=np.ravel(C)   
      
    #calculate erosion for this timeperiod
    E=S*R*K*C #potential erosion
    E[gravel>12.]=E[gravel>12.]-(0.8*E[gravel>12.]) #erosion reduction by 80% if gravel on cropland >12%
    E[E>100.]=100. #max erosion rate is 100t/ha/y
    E = E.reshape(rows,kols) 
     
    #output
    output = Dataset('%s/output/erosion/CC/E_%s_%04i_%s.nc' % (workdir,veget,t,reg),'w')
    output.createDimension('latitude',rows)
    output.createDimension('longitude',kols)
    var=output.createVariable('latitude','f',('latitude',))
    var[:]=[lat0-n*0.0833 for n in range(rows)]
    var=output.createVariable('longitude','f',('longitude',))
    var[:]=[lon0+n*0.0833 for n in range(kols)]
    output.createVariable('E','d',('latitude','longitude',))
    output.variables['E'][:] = E
    output.close()   

if simulation =='LUC':
  ####################################
  # LUC only (CC=set to 1975-1985)
  ###################################
  data0 = Dataset('%s/input_%s/R_factor/ISIMIP2b_R_1975-1985_mean_5m_%s.nc' % (workdir,reg,reg),'r') # R factor
  R = data0.variables['r'][:] 
  R[R<0]=0.;R[R>10e35]=0.
  R[np.isnan(R)==True] = 0.
  R=np.ravel(R)
  data0.close() 

  for t in years:     
    data1 = Dataset('%s/input_%s/C_factor/LUC/C_crop_%04i_5m_%s.nc' % (workdir,reg,t,reg),'r') # C factor
    C = data1.variables['c_lc'][:]
    C[C<0]=0.;C[np.isnan(C)==True] = 0.
    C=np.ravel(C)   
      
    #calculate erosion for this timeperiod
    E=S*R*K*C #potential erosion
    E[gravel>12.]=E[gravel>12.]-(0.8*E[gravel>12.]) #erosion reduction by 80% if gravel on cropland >12%
    E[E>100.]=100. #max erosion rate is 100t/ha/y
    E = E.reshape(rows,kols) 
    
    #output
    output = Dataset('%s/output/erosion/LUC/E_%s_%04i_%s.nc' % (workdir,veget,t,reg),'w')
    output.createDimension('latitude',rows)
    output.createDimension('longitude',kols)
    var=output.createVariable('latitude','f',('latitude',))
    var[:]=[lat0-n*0.0833 for n in range(rows)]
    var=output.createVariable('longitude','f',('longitude',))
    var[:]=[lon0+n*0.0833 for n in range(kols)]
    output.createVariable('E','d',('latitude','longitude',))
    output.variables['E'][:] = E
    output.close()        

