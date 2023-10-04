import tqdm
from specutils.fitting import fit_lines
from specutils.spectra import Spectrum1D
from astropy.modeling import models, fitting
from specutils.fitting import estimate_line_parameters
from astropy.nddata import StdDevUncertainty
from spice_uncertainties import spice_error
from astropy.io import fits

# Choose the key to fit
key = 'Ne VIII 770(Merged)'
cube = data_cube
cube = np.nan_to_num(cube) 

#Computing uncertainties
with fits.open(file) as hdulist:  # specify file name here
    hdu = hdulist[3];               # specify HDU index here
    av_constant_noise_level, sigma = spice_error(hdu, verbose=True)
    sigma = sigma['Total'][0,:,110:750,:].value/10

#Get the array of wavelengths
waves = wavelengths

#Define a function to subtract the background noise
def substract_min_cube(cube):
        det_plane_min = np.nanmin(cube,axis=2)
        for i in range(0,cube.shape[2]): 
            cube[:,:,i] -= det_plane_min
        return cube

#Apply the function
cube_sub = substract_min_cube(cube)


[nx, ny] = cube_sub.shape[1:3]
# Set a reference value for the mean wavelegnth of the emission peak
dopp1= 770

# Initialize arrays
nbrIterConv = np.zeros([nx,ny])
rss, rss1, rss2, rss3 = np.zeros([nx,ny]), np.zeros([nx,ny]), np.zeros([nx,ny]), np.zeros([nx,ny])
amps1, cen1, err1, sig1 = np.zeros([nx,ny]), np.zeros([nx,ny]), np.zeros([nx,ny]), np.zeros([nx,ny])


for i in tqdm.tqdm(range(0,nx)):
    for j in range(0,ny):
        data = np.nan_to_num(cube_sub[:, i, j]).to(u.mW/u.m**2/u.sr/u.Angstrom)*2*2.5 #burn in and spectral binning
        errs = np.nan_to_num(sigma[:, i, j])*1e3
      
        # Create a spectrum object to do the line fitting
        spec = Spectrum1D(flux = data, spectral_axis = waves, uncertainty=StdDevUncertainty(errs))#, mask=mask)
    
        # Create the model : one constant to fit the continuum and one gaussian curve, defining boundaries on the parameters
        c_init = models.Const1D(amplitude = np.nanmean(data)/20, bounds={'amplitude':[np.min(data),np.max(data)]})
        g_init = estimate_line_parameters(spec, models.Gaussian1D(bounds = {'stddev':[0.3,3], 'mean':[dopp1-1, dopp1+1]}))
        
        # Do the line fitting   
        g_fit = fit_lines(spec, g_init +  c_init)
        y_fit = g_fit(waves)

        # Store the results
        amps1[i][j] = (g_fit.amplitude_0.value/2)*np.sqrt(np.pi*2)*g_fit.stddev_0.value
        cen1[i][j] = dopp1  - g_fit.mean_0.value
        rss[i][j] = np.abs(np.mean(g_fit.meta['fit_info']['fvec']))
        nbrIterConv[i][j] = g_fit.meta['fit_info']['nfev']

# Store the result
with open('your_path\\radiance_NeVIII.json', 'wb') as fstd:
    pickle.dump(amps1, fstd)
    
