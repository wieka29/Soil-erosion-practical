# -*- coding: utf-8 -*-
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import *
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import os
import conda

###Do the following only if importing Basemap does not work ###
conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
#check where basemap is installed on your system and modify the following line.
#run in terminal for linux: find `conda info --base` -name epsg
#In my case it's in: '/home/wieka20/anaconda3/share/proj'
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj') 
os.environ["PROJ_LIB"] = proj_lib
###

from mpl_toolkits.basemap import Basemap

# function for moving average
def movingaverage(values,window):
  weights = np.repeat(1.0,window)/window
  smas = np.convolve(values,weights,'valid')
  return smas
#setting
years=range(1975,2006)
workdir='/home/wieka20/Dokumente/Postdoc_LMU/Teaching/Practical/'
reg='BRA' #select your region, 'EU_east' for Eastern Europe, 'EU_west' for Western Europe, 'BRA' for Brazil or 'AUS' for Australia

# input datasets
data= Dataset('%s/input_%s/rest/Global_land_grid_area_5m_grid_ha_%s.nc' % (workdir,reg,reg),'r')
A= data.variables['area'][:] #area of each pixel in ha
rows=len(A[:,0]);kols=len(A[0]) 
lats=np.array(data.variables['lat'][:])
lat0=lats[0];lat1=lats[-1]
lons=np.array(data.variables['lon'][:])
lon0=lons[0];lon1=lons[-1]
llons,llats=np.meshgrid(lons,lats)
A = np.ravel(A)
A[A<0.]=0.
data.close()

######################################
# Exercise 2A: Calculate and plot timeseries
######################################

E_cc_luc=[];E_cc=[];E_luc=[]
for y in range(len(years)):
  #print years[y]
  #erosion CC + LUC
  data = Dataset('%s/output/erosion/E_crop_%04i_%s.nc' % (workdir,years[y],reg),'r')
  E_sum = np.array(data.variables['E'][:])
  E_sum[E_sum==0.]=np.nan
  E_sum[E_sum>100.]=np.nan #exclude unrealistic values
  E_sum=np.ravel(E_sum)*A
  E_cc_luc.append(np.nansum(E_sum)) #global total
  data.close()
  #erosion CC
  data = Dataset('%s/output/erosion/CC/E_crop_%04i_%s.nc' % (workdir,years[y],reg),'r')
  E_sum = np.array(data.variables['E'][:])
  E_sum[E_sum==0.]=np.nan
  E_sum[E_sum>100.]=np.nan 
  E_sum=np.ravel(E_sum)*A
  E_cc.append(np.nansum(E_sum)) 
  data.close()  
  #erosion LUC
  data = Dataset('%s/output/erosion/LUC/E_crop_%04i_%s.nc' % (workdir,years[y],reg),'r')
  E_sum = np.array(data.variables['E'][:])
  E_sum[E_sum==0.]=np.nan
  E_sum[E_sum>100.]=np.nan
  E_sum=np.ravel(E_sum)*A
  E_luc.append(np.nansum(E_sum)) 
  data.close()

E_cc_luc=np.asarray(E_cc_luc)*1E-9 #tonnes to Gt
E_cc=np.asarray(E_cc)*1E-9
E_luc=np.asarray(E_luc)*1E-9

#plotting-moving average
fig, ax = plt.subplots(figsize=(15,10))
fig.subplots_adjust(right=0.75)
ax.plot(years[len(years)-len(movingaverage(E_cc,5)):],movingaverage(E_cc,5),'b-',label='climate change')
ax.plot(years[len(years)-len(movingaverage(E_cc_luc,5)):],movingaverage(E_cc_luc,5),'r-',label='land use change and climate change')
ax.plot(years[len(years)-len(movingaverage(E_luc,5)):],movingaverage(E_luc,5),'g-',label='land use change')
ax.set_xlim([1975,2010])
#ax.set_ylim([2.,5..])
ax.set_xlabel('Calendar year',fontsize=20)
ax.set_ylabel('Soil erosion ($Gt$ $y^{-1}$)',fontsize=20)
ax.set_title('Soil erosion on cropland for %s' % (reg),fontsize=18)
ax.tick_params('x',labelsize=18)
ax.tick_params('y',labelsize=18)
ax.legend(loc='upper left',fontsize=14)
plt.show()
#plt.savefig('%s/output/figure.png' %(workdir))
#plt.close()

######################################
# Exercise 2B:Calculate and plot statistical quantiles for the period 1995-2005 vs 1975-1985
######################################
#calculate,compare and plot quantiles as a boxplot per bin sample data
E0=0.;E1=0.
for y in range(1975,1985):
  print(y)
  #erosion CC + LUC
  data = Dataset('%s/output/erosion/E_crop_%04i_%s.nc' % (workdir,y,reg),'r')
  E = np.array(data.variables['E'][:])
  E[E>20.]=20. #exclude unrealistic values
  E0+=np.ravel(E) #t/ha/y  
  data.close()
for y in range(1995,2005):
  print(y)
  #erosion CC + LUC
  data = Dataset('%s/output/erosion/E_crop_%04i_%s.nc' % (workdir,y,reg),'r')
  E = np.array(data.variables['E'][:])
  E[E>20.]=20.
  E1+=np.ravel(E) 
  data.close()
E_1975_1985=E0/10.;E_1995_2005=E1/10.
E_1975_1985=E_1975_1985[E_1975_1985>0.001];E_1995_2005=E_1995_2005[E_1995_2005>0.001]

data=[E_1975_1985,E_1995_2005]
means=[np.mean(E_1975_1985),np.mean(E_1995_2005)]
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
ax.set_ylabel('Soil erosion ($t$ $ha^{-1}$ $year^{-1}$)',fontsize=20)
ax.set_ylim([0.,5.])
ax.tick_params(axis='y',labelsize=18)
ax.tick_params(axis='x',labelsize=18)
plt.show()
#plt.savefig('%s/output/figure.png' %(workdir)) 
#plt.close()

######################################
#Exercise 2C: plot histograms and compare the distribution of erosion between 1995-2005 and 1975-1985
######################################
fig, ax= plt.subplots()
fig.set_size_inches(10,10)
fig.subplots_adjust(right=0.75)
ax.hist(E_1975_1985, bins=1000, normed=True, alpha=0.5,
         histtype='stepfilled', color='steelblue',
         edgecolor='none')
ax.hist(E_1995_2005, bins=1000, normed=True, alpha=0.5,
         histtype='stepfilled', color='pink',
         edgecolor='none')
ax.set_ylim([0.,1.])
plt.show()
#plt.savefig('%s/output/figure.png' %(workdir)) 
#plt.close()

######################################
#Exercise 3: plot Basemap
######################################
data=(E0/10.).reshape(rows,kols)
fig=plt.figure(figsize=(10,10))
ax=fig.add_subplot(111)
m=Basemap(projection='lcc',resolution='c',lat_0=(lat1-lat0)/2., lon_0=lon0, llcrnrlon=lon0,llcrnrlat=lat1,urcrnrlon=lon1,urcrnrlat=lat0)
m.drawcoastlines()
m.drawparallels(np.arange(lat0,lat1,-10))
m.drawmeridians(np.arange(lon0,lon1,10))
x,y=m(llons,llats)
cmap=plt.cm.inferno_r
data=np.ma.masked_where(data==0.,data)
cmap.set_bad('white')
cs=m.pcolormesh(x,y,data,vmin=0.,vmax=20.,cmap=cmap)
plt.title('Soil erosion 1975-1985 ',fontsize=20)
cbar = m.colorbar(cs,location='bottom',pad="5%",format='%i')
cbar.set_label('$t$ $ha^{-1}$ $year^{-1}$',fontsize=14)
cbar.ax.tick_params(labelsize=14)
plt.show()
#plt.savefig('%s/output/figure.png' %(workdir))  
#plt.close()