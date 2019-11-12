# -*- coding: utf-8 -*-
"""
{

"GW170104":{
      "name":"GW170104",
      "fn_H1"       : "H-H1_LOSC_4_V1-1167559920-32.hdf5",
      "fn_L1"       : "L-L1_LOSC_4_V1-1167559920-32.hdf5",
      "fn_template" : "GW170104_4_template.hdf5",
      "fs"          : 4096,
      "tevent"      : 1167559936.6,
      "utcevent"    : "2017-01-04T10:11:58.60",
      "m1"          : 33.64,
      "m2"          : 24.82,
      "a1"          : -0.236,
      "a2"          : 0.024,
      "approx"      : "lalsim.SEOBNRv2",
      "fband"       : [43.0,800.0],
      "f_min"       : 10.0
  } 

  "GW150914":{
      "name":"GW150914",
      "fn_H1"       : "H-H1_LOSC_4_V2-1126259446-32.hdf5",
      "fn_L1"       : "L-L1_LOSC_4_V2-1126259446-32.hdf5",
      "fn_template" : "GW150914_4_template.hdf5",
      "fs"          : 4096,
      "tevent"      : 1126259462.44,
      "utcevent"    : "2015-09-14T09:50:45.44",
      "m1"          : 41.743,
      "m2"          : 29.237,
      "a1"          : 0.355,
      "a2"          : -0.769,
      "approx"      : "lalsim.SEOBNRv2",
      "fband"       : [43.0,300.0],
      "f_min"       : 10.0
  },
  "LVT151012":{
      "name":"LVT151012",
      "fn_H1"       : "H-H1_LOSC_4_V2-1128678884-32.hdf5",
      "fn_L1"       : "L-L1_LOSC_4_V2-1128678884-32.hdf5",
      "fn_template" : "LVT151012_4_template.hdf5",
      "fs"          : 4096,
      "tevent"      : 1128678900.44,
      "utcevent"    : "2015-10-12T09:54:43.44",
      "m1"          : 44.111,
      "m2"          : 11.205,
      "a1"          : 0.447,
      "a2"          : -0.434,
      "approx"      : "lalsim.SEOBNRv2",
      "fband"       : [43.0,400.0],
      "f_min"       : 10.0
  },
  "GW151226":{
      "name":"GW151226",
      "fn_H1"       : "H-H1_LOSC_4_V2-1135136334-32.hdf5",
      "fn_L1"       : "L-L1_LOSC_4_V2-1135136334-32.hdf5",
      "fn_template" : "GW151226_4_template.hdf5",
      "fs"          : 4096,
      "tevent"      : 1135136350.65,
      "utcevent"    : "2015-12-26T03:38:53.65",
      "m1"          : 19.6427,
      "m2"          : 6.7054,
      "a1"          : 0.3998,
      "a2"          : -0.0396,
      "approx"      : "lalsim.SEOBNRv2",
      "fband"       : [43.0,800.0],
      "f_min"       : 10.0
  },
  
}

"""


import numpy as np 
from matplotlib import pyplot as plt
from matplotlib import gridspec
import h5py
import glob
from scipy.ndimage import gaussian_filter
import os

path="C:/Users/alexk/OneDrive/Documents/phys512/LOSC_Event_tutorial/"
os.chdir(path)

plt.ion()


def read_template(filename):
    dataFile=h5py.File(filename,'r')
    template=dataFile['template']
    th=template[0]
    tl=template[1]
    return th,tl

def read_file(filename):
    dataFile=h5py.File(filename,'r')
    dqInfo = dataFile['quality']['simple']
    qmask=dqInfo['DQmask'][...]

    meta=dataFile['meta']
    gpsStart=meta['GPSstart'].value
    #print meta.keys()
    utc=meta['UTCstart'].value
    duration=meta['Duration'].value
    strain=dataFile['strain']['Strain'].value
    dt=(1.0*duration)/len(strain)

    dataFile.close()
    return strain,dt,utc

def get_fft_vec(n):
    vec=np.arange(n)
    #vec=vec-np.mean(vec)
    vec[vec>n/2]=vec[vec>n/2]-n
    return vec

#def smooth(data0,npix):
#    #make  1-d gaussians of the correct lengths
#    xind=get_fft_vec(data0.size)
#    
#    sig=npix/np.sqrt(8*np.log(2))
#    xvec=np.exp(-0.5*xind**2/sig**2)    
#    #xvec[40000:100000]=0
#    xvec=xvec/xvec.sum()
#    
#    
#    
#    #if we didn't mess up, the kernel FT should be strictly real
#    kernel=np.fft.fftshift(np.real(np.fft.fft(xvec)))
#
#    data_ft=np.fft.fft(data0)
#    data_ft=data_ft*kernel
#    
#    
#    return np.real(np.fft.ifft(data_ft))
    
#arranging the templates in the order that the files are read #mapping
templates=np.append(glob.glob("[G]*.hdf5"),"LVT151012_4_template.hdf5")
template_names=templates.copy()
template_names[0]=templates[2]
template_names[1]=templates[0]
template_names[2]=templates[3]
template_names[3]=templates[1]

#reading the files
fnames=glob.glob("[HL]-*.hdf5")
strain=np.zeros([len(fnames),131072])

for i in range(len(fnames)):
    fname=fnames[i]
    print('reading file ',fname)
    strain[i,:],dt,utc=read_file(fname)


th=np.zeros([len(fnames),131072])
tl=np.zeros([len(fnames),131072])
for i in range(len(template_names)):
    template_name=template_names[i]
    th[i,:],tl[i,:]=read_template(template_name)


x=get_fft_vec(strain[0].size)

#defining the window
window=0.5*(1-np.cos(x*np.pi/np.max(x)))

normfac=np.sqrt(np.mean(window**2))

data=strain*0
dataft=np.fft.rfft(strain)*0
power_spec=dataft*0
data_sm=dataft*0


for i in range(0,strain.shape[0]):
    
    data[i]=(strain[i]*window) #windowing data
    dataft[i]=np.fft.rfft(data[i])
    power_spec[i]=dataft[i]*np.conj(dataft[i])
#Smoothing the data:
    data_sm[i]=(power_spec[i])
    data_sm[i]=gaussian_filter(np.real(data_sm[i]),4)


def get_SNR(t,data,window,data_sm):

    dt=1/4096 #from LIGO
    A=(window*t) #windowing the template
    A_ft=np.fft.rfft(A) #fourier transform of the template
    d=data
    d_ft=np.fft.rfft(d)
    N=data_sm
    N_ift=np.fft.irfft(N)
    #from formula shown in tutorial
    lhalf=np.transpose(N**(-1/2)*A_ft)
    rhalf=(N**(-1/2)*d_ft)

    #finding the match filter
    m=np.fft.fftshift(np.fft.irfft(np.conj(lhalf)*rhalf))
    

    lhalf_ift=np.fft.irfft(lhalf)
    SNR=m*np.sqrt(abs(lhalf_ift)**2) #from tutorial using a matched filter
    SNR_analytic=np.fft.irfft(np.fft.rfft(A)/np.sqrt(data_sm)) #signal divided by noise 
    
    #find the point where half contributes from above and half contributes from below
    power_temp = np.abs(lhalf)**2 
    half_power=np.fft.rfftfreq(len(data),dt)[np.argmin(np.abs(np.cumsum(power_temp)/np.sum(power_temp)-0.5))]
    #finding where the ration of the cumulative sum and the total sum gives 1/2
    #then taking the found frequency at that array position with time step dt 
    return SNR,SNR_analytic,m,half_power




SNR_H=np.zeros((4,data[0].size))
SNR_L=np.zeros((4,data[0].size))
SNR_anal_H=np.zeros((4,data[0].size))
SNR_anal_L=np.zeros((4,data[0].size))
matched_H=np.zeros((4,data[0].size))
matched_L=np.zeros((4,data[0].size))
half_power_H=np.zeros(4)
half_power_L=np.zeros(4)


SNR_plot = plt.figure(figsize=(20,20))
spec = gridspec.GridSpec(ncols=3, nrows=4, figure=SNR_plot, hspace=0.2, wspace=0.3)    
   
#getting data and getting SNR plots
for i in range(0,strain.shape[0]):
    
    if i<4:
        SNR_H[i],SNR_anal_H[i],matched_H[i],half_power_H[i]=get_SNR(th[i],data[i],window,data_sm[i])
        

        SNR_plot.add_subplot(spec[i, 0])
        plt.plot(SNR_H[i])
        plt.title(fnames[i]+' Hanford SNR')
        plt.xlabel("")
        plt.ylabel("SNR_Hanford")
        
        #
        SNR_plot.add_subplot(spec[i, 2])
        plt.plot(((SNR_H[i])**2+(SNR_L[i])**2)**(1/2))
        plt.title(fnames[i]+' Combined SNR')
        plt.xlabel("")
        plt.ylabel("SNR")
        
    else:    
        SNR_L[i-4],SNR_anal_L[i-4],matched_L[i-4],half_power_L[i-4]=get_SNR(tl[i-4],data[i],window,data_sm[i])
        
        SNR_plot.add_subplot(spec[i-4,1])
        plt.plot(SNR_L[i-4])
        plt.title(fnames[i]+' Livingston SNR')
        plt.xlabel("")
        plt.ylabel("SNR_Livingston")
     
plt.savefig("SNR_plots")
matched_plot = plt.figure(figsize=(20,20))   
spec2 = gridspec.GridSpec(ncols=2, nrows=4, figure=matched_plot, hspace=0.2, wspace=0.3)       

#Match filter plots
for i in range(0,strain.shape[0]):        
    if i<4:
        
        matched_plot.add_subplot(spec2[i, 0])
        plt.plot(matched_H[i])
        plt.title(fnames[i]+' Hanford Matched filter')
        plt.xlabel("")
        plt.ylabel("Matched filter")
        
        
    else:         
        
        matched_plot.add_subplot(spec2[i-4, 1])
        plt.plot(matched_L[i-4])
        plt.title(fnames[i]+' Livingston Matched filter')
        plt.xlabel("")
        plt.ylabel("Matched filter")

plt.savefig("Matched_filter")


#for power spectrum plots:

power_plot = plt.figure(figsize=(20,20))   
spec4 = gridspec.GridSpec(ncols=2, nrows=4, figure=power_plot, hspace=0.2, wspace=0.3)   
for i in range(0,strain.shape[0]):        
    if i<4:
        
        power_plot.add_subplot(spec2[i, 0])
        plt.semilogy(data_sm[i])
        plt.title(fnames[i]+' Hanford PS')
        plt.xlabel("")
        plt.ylabel("Power Spectrum")
        
        
    else:         
        
        power_plot.add_subplot(spec2[i-4, 1])
        plt.semilogy(data_sm[i])
        plt.title(fnames[i]+' Livingston PS')
        plt.xlabel("")
        plt.ylabel("Power Spectrum")

plt.savefig("Power_spectrum")

SNR_anal_plot = plt.figure(figsize=(20,20))
spec3 = gridspec.GridSpec(ncols=3, nrows=4, figure=SNR_anal_plot, hspace=0.2, wspace=0.3)    

#analytic SNR plots:
for i in range(0,strain.shape[0]):
    
    if i<4:
       
        SNR_anal_plot.add_subplot(spec3[i, 0])
        plt.plot(SNR_anal_H[i])
        plt.title(fnames[i]+' Hanford analytic SNR')
        plt.xlabel("")
        plt.ylabel("SNR_Hanford")
        
        #
        SNR_anal_plot.add_subplot(spec3[i, 2])
        plt.plot(((SNR_anal_H[i])**2+(SNR_anal_L[i])**2)**(1/2))
        plt.title(fnames[i]+' Combined SNR')
        plt.xlabel("")
        plt.ylabel("SNR")
        
    else:    
       
        
        SNR_anal_plot.add_subplot(spec3[i-4,1])
        plt.plot(SNR_anal_L[i-4])
        plt.title(fnames[i]+' Livingston analytic SNR')
        plt.xlabel("")
        plt.ylabel("SNR_Livingston")

plt.savefig("SNR_analytic")

print("For 1e, the frequencies are for Hanford:"+str(half_power_H ))
print("For 1e, the frequencies are for Livingston:"+str(half_power_L ))

#ANALYSIS FOR F

#finding the time

time=np.arange(len(data[1]))*dt
#for one event find the uncertainty in the time at which the event happened

#after a lot of mental gymnastics, i took the width at 70% (considering the model to be gaussian) and found the uncertainty
matched_unc=np.abs((np.argmax(np.abs(matched_H[0]))-np.argmin(np.abs(np.abs(matched_H[0])-np.max(np.abs(matched_H[0])/1.30))))*dt)

#difference in time between both detectors
difference=np.abs(time[np.argmax(np.abs(SNR_H[0]))]-time[np.argmax(np.abs(SNR_L[0]))])