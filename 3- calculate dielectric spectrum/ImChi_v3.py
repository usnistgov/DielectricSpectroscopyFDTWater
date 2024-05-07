import numpy as np
from multiprocessing import Pool
import matplotlib.pyplot as plt

dipole_file = "dipoleFileReady"    # File with dipoles
t_window = 0.06             # Length of ACF(t) window to use, in ns
t_istart_freq = 0.5E-3      # Separation between time origins to average over, in ns
n_repeat = 20               # Number of resamplings for statistical analysis of autocorrelation version
delta_f = 0.1               # Separation between frequencies for autocorrelation output, in GHz
f_max_autocorr = 1000       # Maximum frequency for autocorrelation, in GHz
f_max_fourier = 5E5         # Maximum frequency for Fourier, in GHz
dlnf_start = 0.1            # Relative frequency resolution for Fourier at low frequency end
dlnf_stop = 0.01            # Relative frequency resolution for Fourier at high frequency end
nf_decade = 100             # Number of frequencies/decade in log part of frequency grid
n_cores = 4                 # Number of cores to parallelize over (single-node only, using multiprocessing)
T = 300                     # K, temperature
#V = 6491.529286             # Average volume of box, in Angstrom^3 for LAMMPS (old SPC/E data)
V = 14768                   # Average volume of box, in Angstrom^3 for LAMMPS (AMOEBA data)
dt = 1.0E-6                 # ns per timestep (consistent with outputs in GHz)

# Units:
q_e = 1.602E-19  # electron charge in coulomb
k_B = 1.380E-23  # Boltzmann constant in J/K
epsilon0 = 8.854e-22  # in Farads/Angstrom

# Load data
P = np.loadtxt(dipole_file) * q_e  # in e*Angstrom
assert P.shape[1] == 3  # no extra columns (else filter/separate out above)

use_Hann_window = False

##########################################

def im_chi_fourier(P: np.ndarray, f: float, df: float) -> tuple[float, float]:
    """
    Compute Im(chi) by Fourier analysis of `P`, at frequency `f`.
    Chunking and accuracy is selected based on frequency resolution `df`.
    Return mean and standard error of Im(chi).
    """
    if use_Hann_window:
        # Hann window for Fourier transform
        t_width = 1.0 / (2 * np.pi * df)  # uncertainty limit
        chunk_spacing = int(np.ceil(t_width/dt))
        t_width = chunk_spacing * dt  # make exact multiple of dt
        chunk_size = 2 * chunk_spacing
        if len(P) <= 2 * chunk_size:
            return np.nan, np.nan
        print(f"Processing frequency {f:.3g}")
        t_chunk = (np.arange(chunk_size) - 0.5*chunk_size) * dt
        kernel = (
            np.exp((-1j * 2 * np.pi * f) * t_chunk) 
            * np.cos(0.5 * np.pi * t_chunk / t_width) ** 2
        )
        t_tot = 0.75 * t_width  # effective time covered by kernel
    else:
        # Gaussian window
        t_sigma = 0.5 / (2 * np.pi * df)  # uncertainty limit
        chunk_spacing = int(np.ceil(2 * t_sigma/dt))
        chunk_size = 3 * chunk_spacing
        if len(P) <= 2 * chunk_size:
            return np.nan, np.nan
        print(f"Processing frequency {f:.3g}")
        t_chunk = (np.arange(chunk_size) - 0.5*chunk_size) * dt
        kernel = np.exp((-1j * 2 * np.pi * f) * t_chunk - 0.5 * (t_chunk / t_sigma)**2)
        t_tot = np.sqrt(np.pi) * t_sigma  # effective time covered by kernel

    # Break data into blocks of chunk_spacing for efficiency:
    n_blocks = len(P) // chunk_spacing
    kernel = kernel.reshape((-1, chunk_spacing))
    P = P[:(n_blocks * chunk_spacing)].reshape((n_blocks, chunk_spacing, -1))
    P_tilde_terms = np.einsum("ac, bci -> abi", kernel, P) * dt

    # Recombine into the Fourier transforms:
    n_kernel_blocks = kernel.shape[0]
    P_tilde = P_tilde_terms[0, : -n_kernel_blocks]
    for i_block in range(1, n_kernel_blocks):
        P_tilde += P_tilde_terms[i_block, i_block : i_block - n_kernel_blocks]
    ImChi = f * (np.abs(P_tilde)**2).sum(axis=-1) * (np.pi / (3 * V * k_B * T * epsilon0 * t_tot))
    return ImChi.mean(), ImChi.std() / np.sqrt(len(ImChi) - 1)


def im_chi_wrapper(args):
    """Wrapper for parallelized execution of im_chi_fourier."""
    return im_chi_fourier(P, *args)

def im_chi_fourier_all(f: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Perform above Fourier analysis over a logarithmic frequency grid `f`.
    The results may not be available for the lowest frquencies,
    and the output frequency grid will be truncated appropriately.
    """
    df = f * np.geomspace(dlnf_start, dlnf_stop, len(f))
    with Pool(n_cores) as pool:
        results = np.array(pool.map(im_chi_wrapper, zip(f, df))).T
    sel = np.where(np.isfinite(results[0]))[0]
    return f[sel], results[0, sel], results[1, sel]


def im_chi_autocorr(P: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute Im(chi) by autocorrelation analysis analysis of `P`.
    Window settings for autocorrelation are specified in global parameters.
    Error analysis is performed by resampling of the contributing windows.
    Return f (GHz), mean Im(chi) and std.err. Im(chi).
    """
    # Downsample time based on needed frequency resolution:
    Ndown = int(np.floor(0.5 / (f_max_autocorr * dt)))  # using f_max as Nyquist freq
    dt_eff = dt * Ndown
    P = P[::Ndown]

    # Time origins for windows:
    n_window = int(t_window/dt_eff)
    n_istart_freq = max(1, int(t_istart_freq/dt_eff))
    n_starts = (len(P) - n_window) // n_istart_freq
    i_starts = np.arange(n_starts) * n_istart_freq

    # Precompute resampling weights for statistical analysis:
    weights = np.random.poisson(1.0, size=(len(i_starts), n_repeat))
    weights = weights / weights.sum(axis=0)[None, :]

    # Compute auto-correlation function averaged (with resampling) over t_start
    autocorr = np.zeros((n_repeat, n_window), dtype=float)
    for i, weight in zip(i_starts, weights):
        autocorr_term = P[i:i+n_window] @ P[i]  # Calculate dot product
        autocorr += np.outer(weight, autocorr_term)  # Average with resampling

    # Scale by denominator to get Eq. 3
    autocorr *= 1.0 / (3 * V * k_B * T * epsilon0)

    # Fourier transform with padding:
    n_FFT_pts = int(np.round(1/(dt_eff*delta_f))) # Convert from delta_f in GHz to n_FFT_pts
    pos_sel = slice(1, n_FFT_pts//2)  # positive frequencies only
    f = np.fft.fftfreq(n_FFT_pts, dt_eff)[pos_sel] # Get corresponding freqs. in 1/ps
    autocorr0 = autocorr[:, :1]
    autocorr_tilde = np.fft.fft(autocorr, n=n_FFT_pts)[:, pos_sel] * dt_eff

    # Switch to chi
    i_omega = 1j * (2*np.pi*f)
    alpha = (autocorr0 - autocorr[:, 1:2]) / (autocorr0 * dt_eff)
    correction_factor = (i_omega + alpha) * dt_eff / (1 - np.exp(-(i_omega + alpha) * dt_eff))
    ImChi = (autocorr0 * correction_factor - i_omega * autocorr_tilde).conj().imag
    return f, ImChi.mean(axis=0), ImChi.std(axis=0)
    

def plot_with_error(x, y_mean, y_std, color, label):
    """Plot 95% (2 sigma) confidence interval."""
    plt.fill_between(x, y_mean - 2*y_std, y_mean + 2*y_std, alpha=0.2, color=color)
    plt.plot(x, y_mean, color=color, label=label)


def combine_results(results_set, f_out):
    """
    Combine multiple result sets weighted by their error.
    Produce output on specified frequency grid.
    """
    yw_sum = np.zeros_like(f_out)
    w_sum = np.zeros_like(f_out)
    for f, y, sigma_y in results_set:
        sel = np.where(np.logical_and(f_out >= f.min(), f_out <= f.max()))[0]
        y_out = np.interp(f_out[sel], f, y)
        w_out = np.interp(f_out[sel], f, sigma_y) ** (-2)
        yw_sum[sel] += y_out * w_out
        w_sum[sel] += w_out
    y_mean = yw_sum / w_sum
    y_std = 1 / np.sqrt(w_sum)
    return f_out, y_mean, y_std


# Set up log frequency grid."""
f_min = delta_f
f_max = min(max(f_max_fourier, f_max_autocorr), np.pi / dt)
f_out = np.geomspace(f_min, f_max, nf_decade * int(np.log10(f_max/f_min)) + 1)

# Get imaginary susceptibility via both methods
results_autocorr = im_chi_autocorr(P)
results_fourier = im_chi_fourier_all(f_out)

# Combine results of both methods
results_combined = combine_results([results_autocorr, results_fourier], f_out)

# Determine Real portion from combined result of imaginary susceptibility
freq = np.array(list(results_combined[:][0]))
ImChi = np.array(list(results_combined[:][1]))
RealChi = np.zeros_like(ImChi)
for i in range(len(ImChi)):
    factor = freq / (freq**2 - freq[i]**2)
    factor[i] = 0
    RealChi[i] = (1/np.pi)*sum((freq[1:]-freq[:-1])*((ImChi[1:]*factor[1:])+(ImChi[:-1]*factor[:-1])))

# Save results
np.savetxt("ImChi.dat", np.array(results_combined[:2]).T)
np.savetxt("RealChi.dat",np.vstack((freq,RealChi)).T)

plt.figure()
plot_with_error(*results_autocorr, "r", "Autocorrelation")
plot_with_error(*results_fourier, "g", "Fourier")
plot_with_error(*results_combined, "k", "Combined")
plt.ylim(1E-3, 100)
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.savefig("ImChi.pdf", bbox_inches="tight")
plt.show()
