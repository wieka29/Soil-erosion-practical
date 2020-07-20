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
#In my case it's in: C:\\Users\\adminuser\\Anaconda2\\Library\\share
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

#set regional settings
if reg=='AUS':
    lat0=-8.1;lat1=-45.1;lon0=111.1; lon1=158.
if reg=='BRA':
    lat0=6.;lat1=-35.;lon0=-75.; lon1=-32.
if reg=='EU_east':
    lat0=82.;lat1=19.85;lon0=11.34;lon1=69.72
if reg=='EU_west':
    lat0=82.;lat1=19.85;lon0=-31.1;lon1=11.34
nc=Dataset('%s/input/Area_ORCHIDEE.nc' %(workdir),'r')
lats=np.array(nc.variables['lat'][:61])
lons=np.array(nc.variables['lon'][:-1])

rows=[]
for i in range(len(lats)):
  if lats[i]>=lat1 and lats[i]<=lat0:
    rows.append(i)
kols=[]
for j in range(len(lons)):
  if lons[j]>=lon0 and lons[j]<=lon1:
    kols.append(j)
rows=np.asarray(rows)
kols=np.asarray(kols)

area=nc.variables['Areas'][rows[0]:rows[-1]+1,kols[0]:kols[-1]+1]
contfrac=nc.variables['CONTFRAC'][rows[0]:rows[-1]+1,kols[0]:kols[-1]+1]
lats_sel=nc.variables['lat'][rows[0]:rows[-1]+1]
lons_sel=nc.variables['lon'][kols[0]:kols[-1]+1]
AREA=area*contfrac
AREA[AREA<0.]=0.;AREA[abs(AREA)>1.e36]=0.
land=contfrac<=1.
llons,llats=np.meshgrid(lons_sel,lats_sel)
nc.close()

######################################
# Exercise 1A: Calculate and plot timeseries
######################################
CE_cc_luc=[];CE_cc=[]
for y in range(len(years)):
  print(years[y])
  #C erosion CC + LUC
  data = Dataset('%s/input/SOC_data/Em_CE_cc_luc_crop_%04i.nc' % (workdir,years[y]),'r')
  CE = np.array(data.variables['E'][rows[0]:rows[-1]+1,kols[0]:kols[-1]+1]) #g/m2
  CE_cc_luc.append(np.nansum(CE*AREA)) #global total
  data.close()
  #C erosion CC
  data = Dataset('%s/input/SOC_data/Em_CE_cc_crop_%04i.nc' % (workdir,years[y]),'r')
  CE = np.array(data.variables['E'][rows[0]:rows[-1]+1,kols[0]:kols[-1]+1]) #g/m2
  CE_cc.append(np.nansum(CE*AREA)) #global total
  data.close()
  
CE_cc_luc_n=np.asarray(CE_cc_luc-CE_cc_luc[0])/np.asarray(CE_cc_luc[0])*100. #normalize
CE_cc_n=np.asarray(CE_cc-CE_cc[0])/np.asarray(CE_cc[0])*100. 
CE_luc=np.asarray(CE_cc[0])+np.asarray(CE_cc_luc)-np.asarray(CE_cc)
CE_luc_n=np.asarray(CE_luc-CE_luc[0])/np.asarray(CE_luc[0])*100. 
#plotting-moving average
fig, ax = plt.subplots(figsize=(15,10))
fig.subplots_adjust(right=0.75)
ax.plot(years[len(years)-len(movingaverage(CE_cc_n,5)):],movingaverage(CE_cc_n,5),'b-',label='climate change')
ax.plot(years[len(years)-len(movingaverage(CE_cc_luc_n,5)):],movingaverage(CE_cc_luc_n,5),'r-',label='land use change and climate change')
ax.plot(years[len(years)-len(movingaverage(CE_luc_n,5)):],movingaverage(CE_luc_n,5),'g-',label='land use change')
ax.set_xlim([1975,2010])
#ax.set_ylim([2.,5..])
ax.set_xlabel('Calendar year',fontsize=20)
ax.set_ylabel('Carbon erosion (% with respect to 1975)',fontsize=20)
ax.set_title('Carbon erosion on cropland for %s' % (reg),fontsize=18)
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
  print(y)
  #erosion CC + LUC
  data = Dataset('%s/input/SOC_data/Em_CE_cc_luc_crop_%04i.nc' % (workdir,y),'r')
  CE = np.array(data.variables['E'][rows[0]:rows[-1]+1,kols[0]:kols[-1]+1]) #g/m2
  CE0+=np.ravel(CE) #g/m2/y  
  data.close()
for y in range(1995,2005):
  print(y)
  #erosion CC + LUC
  data = Dataset('%s/input/SOC_data/Em_CE_cc_luc_crop_%04i.nc' % (workdir,y),'r')
  CE = np.array(data.variables['E'][rows[0]:rows[-1]+1,kols[0]:kols[-1]+1]) #g/m2
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

######################################
#Exercise 2: plot Basemap
######################################
data=(CE0/10.).reshape(len(rows),len(kols))
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
plt.title('Carbon erosion 1975-1985 ',fontsize=20)
cbar = m.colorbar(cs,location='bottom',pad="5%",format='%i')
cbar.set_label('$g$ $m^{-2}$ $year^{-1}$',fontsize=14)
cbar.ax.tick_params(labelsize=14)
plt.show()
#plt.savefig('%s/output/figure.png' %(workdir))  
#plt.close()

