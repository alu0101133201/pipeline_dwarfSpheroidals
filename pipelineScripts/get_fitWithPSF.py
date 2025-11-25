import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys,os
from joblib import Parallel, delayed
from scipy.optimize import curve_fit

def plot_params():
    plt.rcParams['axes.linewidth']    = 2.3
    scatter_kwargs                    = {"zorder":100}
    error_kwargs                      = {"lw":.9, "zorder":0}
    plt.rc('font', weight = 'medium', family = 'serif', size = 15)
    plt.rcParams['xtick.major.size']  = 7
    plt.rcParams['xtick.major.width'] = 1.8
    plt.rcParams['xtick.minor.size']  = 3
    plt.rcParams['xtick.minor.width'] = 1.5
    plt.rcParams['ytick.major.size']  = 7
    plt.rcParams['ytick.major.width'] = 1.8
    plt.rcParams['ytick.minor.size']  = 3
    plt.rcParams['ytick.minor.width'] = 1.8
    plt.rcParams['pdf.fonttype']      = 42
    plt.rcParams['ps.fonttype']       = 42
    plt.rcParams['xtick.direction']   = 'in'
    plt.rcParams['ytick.direction']   = 'in'
    plt.rcParams["legend.fancybox"]   = True
    plt.rcParams['xtick.major.pad']   = 6
    plt.rcParams['ytick.major.pad']   = 6
    plt.rcParams["mathtext.fontset"]  = "dejavuserif"
    return()
plot_params()

star_file=sys.argv[1]
psf_file=sys.argv[2]
star_mag=float(sys.argv[3])
##Test will be made with 1000 values between min and max
back_mean=float(sys.argv[4])
back_std=float(sys.argv[5])
output_folder=sys.argv[6]
ncpu=int(sys.argv[7])
pixScale=float(sys.argv[8])
starSatThresh=float(sys.argv[9])
calFactor=float(sys.argv[10])

def get_profileRange(star_file,star_mag,psf_file,pixScale,starSatTh,calFactor):
    with fits.open(star_file) as hdul:
        I_star=hdul[1].data['SIGCLIP_MEAN']
        R_star=hdul[1].data['RADIUS']
        std_star=hdul[1].data['SIGCLIP_STD']
        mag_star=-2.5*np.log10(I_star*calFactor)+22.5+5*np.log10(pixScale)
        mag_bck=-2.5*np.log10(back_mean*calFactor)+22.5+5*np.log10(pixScale)
        ##First: check where I<6000ADU to avoid saturation and where I>1.5Background
        indexes=np.where((I_star<satThresh)&(~np.isnan(I_star))&(std_star!=0)&(mag_star<mag_bck-0.5))
        indexes_belBck=np.where((~np.isnan(I_star))&(std_star!=0)&(mag_star>mag_bck+1.0))
        if len(indexes[0])==0:
            #If no values are found, probably the first values taken are close to the background, so we
            #make a range less restrictive
            indexes=np.where((I_star<satThresh)&(~np.isnan(I_star))&(std_star!=0)&(mag_star<mag_bck-0.2))
            
        if (len(indexes[0])==0):
            #If still no values are found, then we can assume scatter light won't affect the image
            #So, we give Rmin=100 and Rmax=50 in order to return a scale of 0 on the later on sanity check
            Rmin=100
            Rmax=50
            indexes_ok=[[0]]
        else:
            #We take the first and last value of the indexes
            Rmin=R_star[indexes[0][0]]
            Rmax=R_star[indexes_belBck[0][0]] if len(indexes_belBck[0])>0 else R_star[indexes[0][-1]]
            indexes_ok=indexes[0][(R_star[indexes[0]]>=Rmin)&(R_star[indexes[0]]<=Rmax)]
    return Rmin,Rmax,indexes_ok 

#r_min,r_max=get_profileRange(star_file,star_mag,psf_file,pixScale)




def get_parameterRange(star_file,psf_file,r_min,r_max,back_mean):
    ###BACKGROUND: an uncertainty of 10%
    back_min=0.99*back_mean
    back_max=1.01*back_mean
    ####ALFA: first we compute a very rough value
    table_star=fits.open(star_file)[1].data
    table_psf=fits.open(psf_file)[1].data
    R_star_all=table_star['RADIUS']
    I_star_all=table_star['SIGCLIP_MEAN']
    I_psf_all=table_psf['MEAN']
    indexes=np.where((R_star_all>=r_min)&(R_star_all<=r_max)&(~np.isnan(I_star_all)))
    alf_mean=np.mean((I_star_all[indexes]-back_mean)/I_psf_all[indexes])
    #Second: we give an uncertainty of 20%
    alf_min=0.7*alf_mean
    alf_max=1.3*alf_mean
    if alf_mean<200 and alf_mean>1:
        alf_min=1
        alf_max=200
    
    return(back_min,back_max,alf_min,alf_max,alf_mean)

#back_min,back_max,alfa_min,alfa_max,alfa_mean=get_parameterRange(star_file,psf_file,r_min,r_max,back_mean)

def chisqi(Istar,stdstar,Ipsf,alfa,back):
    """
    Function to compute the i value of the sum inside X2
    """
    Itheo=alfa*Ipsf+back
    num=Istar-Itheo
    return((num**2)/(stdstar**2))

def chisq(star_file,psf_file,alfa,back,r_min,r_max):
    table_star=fits.open(star_file)[1].data
    table_psf=fits.open(psf_file)[1].data
    R_star_all=table_star['RADIUS']
    I_star_all=table_star['SIGCLIP_MEAN']
    std_star_all=table_star['SIGCLIP_STD']
    I_psf_all=table_psf['MEAN']
    ##We assume that R_star==R_psf
    indexes=np.where((R_star_all>=r_min)&(R_star_all<=r_max))
    chi=0
    for index in indexes[0]:
        if np.isnan(I_star_all[index]):
            chi_i=0
        else:
            chi_i=chisqi(I_star_all[index],std_star_all[index],I_psf_all[index],alfa,back)
        chi+=chi_i
    return chi

def compute_chi2_parallel_joblib(star_file,psf_file,alfa_min,alfa_max,back_min,back_max,r_min,r_max,n_points=200):
    alfa_test = np.linspace(alfa_min, alfa_max, n_points)
    back_test = np.linspace(back_min, back_max, n_points)
    chi2test = np.zeros((n_points, n_points))
    def compute_single_chi2(i, j):
        return (j, i, chisq(star_file, psf_file, alfa_test[i], back_test[j], r_min, r_max))
    indices = [(i, j) for i in range(alfa_test.size) for j in range(back_test.size)]
    
    # Use joblib to parallelize
    results = Parallel(n_jobs=ncpu)(delayed(compute_single_chi2)(i, j) for i, j in indices)
    
    # Fill the results array
    for j, i, chi_val in results:
        chi2test[j, i] = chi_val
        
    return chi2test, alfa_test, back_test

def linear_model(I_psf,alfa,back):
    return alfa*I_psf+back


#alfa_test=np.linspace(alfa_min,alfa_max,200)
#back_test=np.linspace(back_min,back_max,200)
#chi2test=np.zeros((200,200))
#for i in range(alfa_test.size):
#    for j in range(back_test.size):
#        chis=chisq(star_file,psf_file,alfa_test[i],back_test[j],r_min,r_max)
#        chi2test[j,i]=chis
#chi2test,alfa_test,back_test=compute_chi2_parallel_joblib(star_file,psf_file,alfa_min,alfa_max,back_min,back_max,r_min,r_max)
#chi2test_min=np.amin(chi2test[~np.isnan(chi2test)])
#indices_chi2test=np.where(chi2test==chi2test_min)
#alfa_min_test=alfa_test[indices_chi2test[1]]
#back_min_test=back_test[indices_chi2test[0]]
#
#aa,bb=np.meshgrid(alfa_test,back_test)
#levels=[0,2.3,6.18,11.83]
#dchi=chi2test-chi2test_min
#base=os.path.basename(star_file)
#frameNumber=base.split('.')[0].split('_')[2]
#plt.figure(figsize=(10,8))
#plt.contour(aa,bb,dchi,levels,colors='k')
#plt.contourf(aa,bb,dchi,levels)
#plt.plot(alfa_min_test,back_min_test,color='r',marker='X',markersize=5,linewidth=0,label=rf'Min:$\alpha$={alfa_min_test}; Bck.={back_min_test}')
#plt.xlabel(r'$\alpha$')
#plt.ylabel('Background')
#plt.title(fr"$\Delta\chi^{2}$: Frame {frameNumber} R: {r_min} - {r_max} ")
#plt.legend()
#plt.tight_layout()
#plt.savefig(f'{output_folder}/{base}_fit.pdf')
#
#print(alfa_min_test[0])

##We are gonna use scipy curve fitting for make the fit

starTab=fits.open(star_file)[1].data
psfTab=fits.open(psf_file)[1].data
R_star=starTab['RADIUS']
I_star=starTab['SIGCLIP_MEAN']
std_star=starTab['SIGCLIP_STD']
###Sanity check: for stars where the first nonan values are compatible with pure background, we will directly return a 0 value for the scale
first_nonans=I_star[~np.isnan(I_star)][:20]

if (np.mean(first_nonans)<=1.005*back_mean)and(np.mean(first_nonans)>=0.995*back_mean):
    alfa_fit=0.000
else:
    r_min,r_max,indexes=get_profileRange(star_file,back_mean,back_std,starSatTh,calFactor,pixScale)
    
    if r_min>=r_max:
        #If for some resaons R_min=R_max (i.e, only one value escape our conditions), we can assume is just an statistics anomally, so we set the scale to 0
        alfa_fit=0.000
    else:
        #NOTE: R_psf might star in 0 and R_star in R_psf[1]. We need to check
        if R_star[0]>psfTab['RADIUS'][0]:
            #We need to remove the first value of the psf
            psfTab=psfTab[1:]
        
        I_psf_ok=psfTab['MEAN'][indexes]
        I_star_ok=I_star[indexes]
        std_star_ok=std_star[indexes]
        
        lstd_star_ok=std_star_ok/(np.log(10)*I_star_ok)
        back_min,back_max,alfa_min,alfa_max,alfa_mean=get_parameterRange(star_file,psf_file,r_min,r_max,back_mean)
        initial_guess=[alfa_mean,back_mean]
        bounds=([alfa_min,back_min],[alfa_max,back_max])
        #Last sanity check: if background is dominant, the scale will be an strange value, probably negative. 
        # To avoid that, we will again give a scale of 0 when this happens
        if (alfa_mean<0) or (alfa_min>alfa_max):
            alfa_fit=0.000
        else:
        
            popt,pcov=curve_fit(linear_model,I_psf_ok,np.log10(I_star_ok),p0=initial_guess,bounds=bounds,sigma=lstd_star_ok,absolute_sigma=True)
            alfa_fit,back_fit=popt
            perr=np.sqrt(np.diag(pcov))

print(alfa_fit)