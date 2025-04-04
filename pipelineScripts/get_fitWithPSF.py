import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys,os
from joblib import Parallel, delayed

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
r_min=float(sys.argv[3])
r_max=float(sys.argv[4])
star_mag=float(sys.argv[5])
##Test will be made with 1000 values between min and max
back_mean=float(sys.argv[6])
output_folder=sys.argv[7]
ncpu=int(sys.argv[8])

back_min=back_mean-500
back_max=back_mean+500

#THIS NUMBERS ARE FOR LBT G PSF, WHOSE MEAN IS 77; BE AWARE OF THAT
if star_mag>=7 and star_mag<8:
    alfa_min=10000
    alfa_max=20000
elif star_mag>=8 and star_mag<9:
    alfa_min=6000 
    alfa_max=12000
elif star_mag>=9 and star_mag<10:
    alfa_min=2000
    alfa_max=8000
elif star_mag>=10 and star_mag<11:
    alfa_min=500
    alfa_max=3000
elif star_mag


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

#alfa_test=np.linspace(alfa_min,alfa_max,200)
#back_test=np.linspace(back_min,back_max,200)
#chi2test=np.zeros((200,200))
#for i in range(alfa_test.size):
#    for j in range(back_test.size):
#        chis=chisq(star_file,psf_file,alfa_test[i],back_test[j],r_min,r_max)
#        chi2test[j,i]=chis
chi2test,alfa_test,back_test=compute_chi2_parallel_joblib(star_file,psf_file,alfa_min,alfa_max,back_min,back_max,r_min,r_max)
chi2test_min=np.amin(chi2test[~np.isnan(chi2test)])
indices_chi2test=np.where(chi2test==chi2test_min)
alfa_min_test=alfa_test[indices_chi2test[1]]
back_min_test=back_test[indices_chi2test[0]]

aa,bb=np.meshgrid(alfa_test,back_test)
levels=[0,2.3,6.18,11.83]
dchi=chi2test-chi2test_min
base=os.path.basename(star_file)
frameNumber=base.split('_')[2]
hNumber=base.split('_')[-1].split('.')[0][1]
plt.figure(figsize=(10,8))
plt.contour(aa,bb,dchi,levels,colors='k')
plt.contourf(aa,bb,dchi,levels)
plt.plot(alfa_min_test,back_min_test,color='r',marker='X',markersize=5,linewidth=0,label=rf'Min:$\alpha$={alfa_min_test}; Bck.={back_min_test}')
plt.xlabel(r'$\alpha$')
plt.ylabel('Background')
plt.title(fr"$\Delta\chi^{2}$: Frame {frameNumber}, CCD{hNumber}")
plt.legend()
plt.tight_layout()
plt.savefig(f'{output_folder}/{base}_fit.pdf')
