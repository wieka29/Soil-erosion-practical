# -*- coding: utf-8 -*-
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import *
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import os
import conda
conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
#check where basemap is installed on your system and modify the following line.
#In my case it's in: C:\\Users\\adminuser\\Anaconda2\\Library\\share
proj_lib = os.path.join(os.path.join(conda_dir, 'Library'), 'share') 
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
from PythonTools import *

# function for moving average
def movingaverage(values,window):
  weights = np.repeat(1.0,window)/window
  smas = np.convolve(values,weights,'valid')
  return smas
#setting
years=range(1975,2006)
workdir='C:\VIKA\Postdoc_LSCE\Indices'
clim_class=("A","B","C","D")
clim_sel=1 #zone A
rows=61;kols=96

data = Dataset('%s\input\climate_classification_aggr_coarse.nc' % (workdir) ,'r') #climate zone data
clim= data.variables['climate'][:61,:] #1= Zone A (tropical), 2= Zone B (arid/dry), 3= Zone C (temperate), 4= Zone D (Cold)
clim[clim<0.]=np.nan
data.close()  

nc=Dataset('%s/input/AREA_ORCHIDEE.nc' %(workdir),'r')
area=nc.variables['Areas'][:61,:-1]
contfrac=nc.variables['CONTFRAC'][:61,:-1]
AREA=area*contfrac
AREA[AREA<0.]=0.;AREA[abs(AREA)>1.e36]=0.
land=contfrac<=1.
lats=np.array(nc.variables['lat'][:61])
lons=np.array(nc.variables['lon'][:-1])
llons,llats=np.meshgrid(lons,lats)
nc.close()

######################################
# Exercise 1A: Calculate and plot timeseries
######################################
CE_cc_luc=[];CE_cc=[]
for y in range(len(years)):
  print years[y]
  #C erosion CC + LUC
  data = Dataset('%s/input/SOC_data/Em_CE_cc_luc_crop_%04i.nc' % (workdir,years[y]),'r')
  CE = np.array(data.variables['E'][:]) #g/m2
  CE[clim!=clim_sel]=0.
  CE_cc_luc.append(np.nansum(CE*AREA)) #global total
  data.close()
  #C erosion CC
  data = Dataset('%s/input/SOC_data/Em_CE_cc_crop_%04i.nc' % (workdir,years[y]),'r')
  CE = np.array(data.variables['E'][:]) #g/m2
  CE[clim!=clim_sel]=0.
  CE_cc.append(np.nansum(CE*AREA)) #global total
  data.close()
  
CE_cc_luc=np.asarray(CE_cc_luc)*1E-15 #g to Pg
CE_cc=np.asarray(CE_cc)*1E-15
CE_luc=CE_cc[0]+(CE_cc_luc-CE_cc)

#plotting-moving average
fig, ax = plt.subplots(figsize=(15,10))
fig.subplots_adjust(right=0.75)
ax.plot(years[len(years)-len(movingaverage(CE_cc,5)):],movingaverage(CE_cc,5),'b-',label='climate change')
ax.plot(years[len(years)-len(movingaverage(CE_cc_luc,5)):],movingaverage(CE_cc_luc,5),'r-',label='land use change and climate change')
ax.plot(years[len(years)-len(movingaverage(CE_luc,5)):],movingaverage(CE_luc,5),'g-',label='land use change')
ax.set_xlim([1975,2010])
#ax.set_ylim([2.,5..])
ax.set_xlabel('Calendar year',fontsize=20)
ax.set_ylabel('Carbon erosion ($Gt C$ $y^{-1}$)',fontsize=20)
ax.set_title('Carbon erosion on cropland for zone %s' % (clim_class[clim_sel-1]),fontsize=18)
ax.tick_params('x',labelsize=18)
ax.tick_params('y',labelsize=18)
ax.legend(loc='middle right',fontsize=14)
plt.show()
#plt.savefig('%s/output/figure.png' %(workdir))
#plt.close()

######################################
# Exercise 1B:Calculate and plot statistical quantiles for the period 1995-2005 vs 1975-1985
######################################
#calculate,compare and plot quantiles as a boxplot per bin sample data
CE0=0.;CE1=0.
for y in range(1975,1985):
  print y
  #erosion CC + LUC
  data = Dataset('%s/input/SOC_data/Em_CE_cc_luc_crop_%04i.nc' % (workdir,y),'r')
  CE = np.array(data.variables['E'][:]) #g/m2
  CE[clim!=clim_sel]=0.  
  CE0+=np.ravel(CE) #g/m2/y  
  data.close()
for y in range(1995,2005):
  print y
  #erosion CC + LUC
  data = Dataset('%s/input/SOC_data/Em_CE_cc_luc_crop_%04i.nc' % (workdir,y),'r')
  CE = np.array(data.variables['E'][:]) #g/m2
  CE[clim!=clim_sel]=0.  
  CE1+=np.ravel(CE) 
  data.close()
CE_1975_1985=CE0/10.;CE_1995_2005=CE1/10.
CE_1975_1985=CE_1975_1985[CE_1975_1985>0.0001];CE_1995_2005=CE_1995_2005[CE_1995_2005>0.0001]

data=[CE_1975_1985,CE_1995_2005]
means=[np.mean(CE_1975_1985),np.mean(CE_1995_2005)]
fig, ax= plt.subplots()
fig.set_size_inches(10,10)
fig.subplots_adjust(right=0.75)
ax.yaxis.grid(True,linestyle='-',which='major',color='lightgrey',alpha=0.5)
ax.set_axisbelow(True)
bp=ax.boxplot(data,notch=0,vert=1,whis=2.5)
#ax.set_yscale('log')
for i in range(2):
  if i==0:
    plt.setp(bp['boxes'][i],color='g')
    plt.setp(bp['whiskers'][i],color='g')
    #plt.setp(bp['fliers'][i],color='g',marker='+')
  else:
    plt.setp(bp['boxes'][i],color='k')
    plt.setp(bp['whiskers'][i],color='k')
    #plt.setp(bp['fliers'][i],color='k',marker='+')
for median in bp['medians']:
  median.set(color='r',linewidth=2)
for whisker in bp['whiskers']:
  whisker.set(color='k',linestyle='--')
xtickNames=plt.setp(ax,xticklabels=np.array(['1975-1985','1995-2005']))
means_all=[means[0],means[1]]
plt.scatter(np.arange(1,3),means_all,color='b',marker='*',s=60)
ax.set_ylabel('Carbon erosion ($g$ $m^{-2}$ $year^{-1}$)',fontsize=20)
ax.set_ylim([0.,7.5])
ax.tick_params(axis='y',labelsize=18)
ax.tick_params(axis='x',labelsize=18)
plt.show()
#plt.savefig('%s/output/figure.png' %(workdir)) 
#plt.close()

#######################################
##Exercise 1C: plot histograms and compare the distribution of erosion between 1995-2005 and 1975-1985
#######################################
#fig, ax= plt.subplots()
#fig.set_size_inches(10,10)
#fig.subplots_adjust(right=0.75)
#ax.hist(CE_1975_1985, bins=100, normed=True, alpha=0.5,
#         histtype='stepfilled', color='steelblue',
#         edgecolor='none')
#ax.hist(CE_1995_2005, bins=100, normed=True, alpha=0.5,
#         histtype='stepfilled', color='pink',
#         edgecolor='none')
#ax.set_ylim([0.,1.])
#plt.show()
##plt.savefig('%s/output/figure.png' %(workdir)) 
##plt.close()

######################################
#Exercise 2: plot Basemap
######################################
data=(CE0/10.).reshape(rows,kols)
fig=plt.figure(figsize=(10,10))
ax=fig.add_subplot(111)
m=Basemap(projection='moll',resolution='c', lon_0=0, llcrnrlon=-180.,llcrnrlat=-60.,urcrnrlon=180.,urcrnrlat=75.)
m.drawcoastlines()
m.drawparallels(np.arange(-60.,75.,30))
m.drawmeridians(np.arange(-180.,180.,30))
x,y=m(llons,llats)
cmap=plt.cm.inferno_r
data=np.ma.masked_where(data==0.,data)
cmap.set_bad('white')
cs=m.pcolormesh(x,y,data,vmin=0.,vmax=20.,cmap=cmap)
plt.title('Carbon erosion 1975-1985 ',fontsize=20)
cbar = m.colorbar(cs,location='bottom',pad="5%",format='%i')
cbar.set_label('$g$ $m^{-2}$ $year^{-1}$',fontsize=14)
cbar.ax.tick_params(labelsize=14)
plt.show()
#plt.savefig('%s/output/figure.png' %(workdir))  
#plt.close()

######################################
#Exercise 3: relationship between soil and C erosion rates 
######################################
CE0=0.;CE1=0.
for y in range(1975,1985):
  #print y
  #erosion CC + LUC
  data = Dataset('%s/input/SOC_data/Em_CE_cc_luc_crop_%04i.nc' % (workdir,y),'r')
  CE = np.array(data.variables['E'][:]) #g/m2
  CE[clim!=clim_sel]=0.  
  CE0+=np.ravel(CE) #g/m2/y  
  data.close()
for y in range(1995,2005):
  #print y
  #erosion CC + LUC
  data = Dataset('%s/input/SOC_data/Em_CE_cc_luc_crop_%04i.nc' % (workdir,y),'r')
  CE = np.array(data.variables['E'][:]) #g/m2
  CE[clim!=clim_sel]=0.  
  CE1+=np.ravel(CE) 
  data.close()

E0=0.;E1=0.
for y in range(1975,1985):
  #print y
  #erosion CC + LUC
  data = Dataset('%s/output/erosion/E_crop_%04i_zone%s.nc' % (workdir,y,clim_class[clim_sel-1]),'r')
  E = np.array(data.variables['E'][:,:-2])
  E[E>20.]=20. #exclude unrealistic values
  E0+=E #t/ha/y  
  data.close()
for y in range(1995,2005):
  #print y
  #erosion CC + LUC
  data = Dataset('%s/output/erosion/E_crop_%04i_zone%s.nc' % (workdir,y,clim_class[clim_sel-1]),'r')
  E = np.array(data.variables['E'][:,:-2])
  E[E>20.]=20.
  E1+=E
  data.close()
E_1975_1985=E0/10.;E_1995_2005=E1/10.

data0=1.*E_1975_1985
data0[np.isnan(data0)==True]=0.

data0=np.ma.vstack((np.ma.masked_all((48,4320)),data0,np.ma.masked_all((432,4320)))) # fill pixels (0.083 in lat) north of 83.9917 and pixels south of -55.9523
data0_aggr=regrid(data0[:],[90.,-90.],1./12.,[90.,-90.],2.5,1/12.,[-180.,180.],1./12.,[-180.,180.],3.75,1/24.,'d','d',True,True,True,1)

E0_data=1*np.ravel(data0_aggr[:61,:-1])
CE0_data=(CE0/10.)
E0_data[E0_data==0.]=np.nan
CE0_data[CE0_data==0.]=np.nan
E0_data[np.isnan(CE0_data)==True]=np.nan
CE0_data[np.isnan(E0_data)==True]=np.nan
E0_data=E0_data[np.isnan(E0_data)==False]
CE0_data=CE0_data[np.isnan(CE0_data)==False]

data=1.*E_1995_2005
data[np.isnan(data)==True]=0.

data=np.ma.vstack((np.ma.masked_all((48,4320)),data,np.ma.masked_all((432,4320)))) # fill pixels (0.083 in lat) north of 83.9917 and pixels south of -55.9523
data_aggr=regrid(data[:],[90.,-90.],1./12.,[90.,-90.],2.5,1/12.,[-180.,180.],1./12.,[-180.,180.],3.75,1/24.,'d','d',True,True,True,1)

E_data=1*np.ravel(data_aggr[:61,:-1])
CE_data=(CE1/10.)
E_data[E_data==0.]=np.nan
CE_data[CE_data==0.]=np.nan
E_data[np.isnan(CE_data)==True]=np.nan
CE_data[np.isnan(E_data)==True]=np.nan
E_data=E_data[np.isnan(E_data)==False]
CE_data=CE_data[np.isnan(CE_data)==False]

fig, ax = plt.subplots()
fig.subplots_adjust(right=0.75)
z0 = np.polyfit(E0_data,CE0_data, 1)
z = np.polyfit(E_data,CE_data, 1)
p0=np.poly1d(z0)
p=np.poly1d(z)
ax.plot(E0_data, p(E0_data), 'k-')
ax.plot(E_data, p(E_data), 'r-')
ax.scatter(E0_data,CE0_data,color='k')
ax.scatter(E_data,CE_data,color='r')
ax.set_xlabel('Soil erosion ($t$ $ha^{-1}$ $year^{-1}$)',fontsize=22)
ax.set_ylabel('C erosion ($g$ $m^{-2}$ $year^{-1}$)',fontsize=22)
ax.tick_params(axis='y',labelsize=22)
ax.tick_params(axis='x',labelsize=22)
#ax.legend(loc='upper left',fontsize=20)
plt.show()
