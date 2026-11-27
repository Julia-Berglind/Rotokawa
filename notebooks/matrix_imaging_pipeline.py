"""
matrix_imaging_pipeline.py
First-pass passive matrix imaging from ambient seismic noise (vertical component).

Pipeline
  1. load          ObsPy Stream -> array [n_sta, n_t], station coords in local km
  2. preprocess    detrend, bandpass, running-absolute-mean temporal normalisation
  3. spectra       cut into overlapping windows, FFT, spectral whitening -> U[w, i, f]
  4. K(f)          response matrix  K_ij(f) = < U_i(f) U_j(f)* >  (5x5 per frequency)
                   its inverse FFT gives the noise cross-correlations C_ij(tau)
  5. gate          mute the direct inter-station surface wave (keep later "scattered" lags)
  6. focus         project K onto a grid of virtual sources/receivers with an assumed velocity:
                     R(r_in, r_out) = sum_f sum_ij exp(+i w d_i(r_out)/v) K_ij exp(+i w d_j(r_in)/v)
                   confocal image I(r) = |R(r, r)|, plus the focused matrix around one point
                   and the array point-spread function (PSF) for a synthetic point scatterer.

Usage (in the environment that has ObsPy):
    from obspy import read
    import matrix_imaging_pipeline as mip
    st = read("NZ.*.EHZ.sac")          # or the processed Stream already in memory
    out = mip.run(st)                  # writes figures + mi_results.npz to the current folder

Parameters that change the result (all defaults are my picks, not established values):
    fmin, fmax   correlation/imaging band (Hz)            default 0.1-1.0
    v            assumed surface-wave velocity (km/s)      default 2.0
    win_s        correlation window length (s)            default 120
    v_mute       direct-wave mute velocity (km/s)          default 1.0  (mutes |tau| < d/v_mute + t_pad)
"""
import numpy as np
from scipy.signal import butter, sosfiltfilt, detrend
from scipy.ndimage import uniform_filter1d

KM_PER_DEG = 111.195


# ---------------------------------------------------------------- 1. load
def stream_to_array(st, coords=None):
    """coords: optional {station: (lat, lon)}; otherwise read from tr.stats.sac.stla/stlo."""
    st = st.copy()
    st.sort(keys=["station"])
    names = [tr.stats.station for tr in st]
    fs = float(st[0].stats.sampling_rate)
    if any(tr.stats.sampling_rate != fs for tr in st):
        raise ValueError("all traces must share one sampling rate")
    n = min(tr.stats.npts for tr in st)
    data = np.array([tr.data[:n].astype(float) for tr in st])
    if coords is None:
        lat = np.array([tr.stats.sac.stla for tr in st])
        lon = np.array([tr.stats.sac.stlo for tr in st])
    else:
        lat = np.array([coords[s][0] for s in names])
        lon = np.array([coords[s][1] for s in names])
    return names, fs, data, lat, lon


def to_local_km(lat, lon, lat0=None, lon0=None):
    lat0 = np.mean(lat) if lat0 is None else lat0
    lon0 = np.mean(lon) if lon0 is None else lon0
    x = (lon - lon0) * np.cos(np.radians(lat0)) * KM_PER_DEG
    y = (lat - lat0) * KM_PER_DEG
    return x, y, (lat0, lon0)


# ---------------------------------------------------------------- 2. preprocess
def preprocess(data, fs, fmin=0.1, fmax=1.0, norm_win_s=None):
    """Detrend, bandpass (4th-order Butterworth, zero phase), running-absolute-mean normalisation
    (Bensen et al. 2007 style). norm_win_s defaults to half the longest period."""
    d = detrend(data, axis=1, type="linear")
    sos = butter(4, [fmin, fmax], btype="bandpass", fs=fs, output="sos")
    d = sosfiltfilt(sos, d, axis=1)
    if norm_win_s is None:
        norm_win_s = 0.5 / fmin
    nw = max(1, int(round(norm_win_s * fs)))
    env = uniform_filter1d(np.abs(d), size=nw, axis=1)
    return d / (env + 1e-12 * env.max())


# ---------------------------------------------------------------- 3. windowed whitened spectra
def band_taper(f, fmin, fmax, frac=0.1):
    """Cosine-tapered box: 0 outside [fmin, fmax], 1 inside, tapers over frac of the band at each edge."""
    w = np.zeros_like(f)
    bw = (fmax - fmin) * frac
    inside = (f >= fmin) & (f <= fmax)
    w[inside] = 1.0
    lo = (f >= fmin) & (f < fmin + bw)
    hi = (f > fmax - bw) & (f <= fmax)
    w[lo] = 0.5 * (1 - np.cos(np.pi * (f[lo] - fmin) / bw))
    w[hi] = 0.5 * (1 - np.cos(np.pi * (fmax - f[hi]) / bw))
    return w


def window_spectra(data, fs, win_s=120.0, overlap=0.5, fmin=0.1, fmax=1.0, smooth_hz=0.02):
    """Returns U[w, i, f] (whitened, band-limited), f, nfft. Zero-padded to 2x window to avoid wrap-around."""
    L = int(round(win_s * fs))
    step = max(1, int(round(L * (1 - overlap))))
    starts = np.arange(0, data.shape[1] - L + 1, step)
    nfft = 2 * L
    f = np.fft.rfftfreq(nfft, 1 / fs)
    df = f[1] - f[0]
    taper = np.hanning(L)
    bt = band_taper(f, fmin, fmax)
    nsm = max(1, int(round(smooth_hz / df)))
    U = np.empty((len(starts), data.shape[0], len(f)), complex)
    for k, s in enumerate(starts):
        seg = detrend(data[:, s:s + L], axis=1) * taper
        S = np.fft.rfft(seg, n=nfft, axis=1)
        amp = uniform_filter1d(np.abs(S), size=nsm, axis=1)
        U[k] = S / (amp + 1e-12 * amp.max()) * bt
    return U, f, nfft


# ---------------------------------------------------------------- 4. response matrix and correlations
def response_matrix(U):
    """K[f, i, j] = mean over windows of U_i U_j^*  (Hermitian 5x5 at each frequency)."""
    return np.einsum("wif,wjf->fij", U, U.conj()) / U.shape[0]


def correlations(K, nfft, fs):
    """C[i, j, tau] = sum_t u_i(t+tau) u_j(t); positive lag = wave travelling j -> i."""
    C = np.fft.irfft(np.transpose(K, (1, 2, 0)), n=nfft, axis=-1)
    C = np.fft.fftshift(C, axes=-1)
    lags = (np.arange(nfft) - nfft // 2) / fs
    return C, lags


def correlations_to_K(C, nfft):
    Cs = np.fft.ifftshift(C, axes=-1)
    return np.transpose(np.fft.rfft(Cs, n=nfft, axis=-1), (2, 0, 1))


# ---------------------------------------------------------------- 5. gate the direct wave
def distance_matrix(x, y):
    return np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])


def gate(C, lags, dist, v_mute=1.0, t_pad=3.0, t_max=60.0, edge_s=2.0):
    """Zero |tau| < dist/v_mute + t_pad (direct wave, and the zero-lag spike on the diagonal)
    and |tau| > t_max, with cosine edges."""
    a = np.abs(lags)
    G = np.empty_like(C)
    for i in range(C.shape[0]):
        for j in range(C.shape[1]):
            t0 = dist[i, j] / v_mute + t_pad
            w = np.clip((a - t0) / edge_s, 0, 1)
            w = 0.5 * (1 - np.cos(np.pi * w))
            w *= np.clip((t_max - a) / edge_s, 0, 1)
            G[i, j] = C[i, j] * w
    return G


# ---------------------------------------------------------------- 6. focusing
def make_grid(x, y, pad_km=10.0, dx_km=0.5):
    gx = np.arange(x.min() - pad_km, x.max() + pad_km + dx_km / 2, dx_km)
    gy = np.arange(y.min() - pad_km, y.max() + pad_km + dx_km / 2, dx_km)
    GX, GY = np.meshgrid(gx, gy)
    return gx, gy, GX.ravel(), GY.ravel()


def steering(f, x, y, px, py, v):
    """E[f, i, p] = exp(+i 2 pi f d_ip / v): undoes the travel-time phase from station i to pixel p."""
    d = np.hypot(x[:, None] - px[None, :], y[:, None] - py[None, :])  # [i, p]
    return np.exp(2j * np.pi * f[:, None, None] * d[None] / v)


def confocal_image(K, E, fmask):
    """I(p) = | sum_f sum_ij E_ip K_ij E_jp |  (virtual source = virtual receiver = p)."""
    Kf, Ef = K[fmask], E[fmask]
    return np.abs(np.einsum("fip,fij,fjp->p", Ef, Kf, Ef))


def focused_row(K, E, fmask, p_in):
    """|R(r_in = p_in, r_out = every pixel)|: where energy focused at p_in comes back to."""
    Kf, Ef = K[fmask], E[fmask]
    return np.abs(np.einsum("fip,fij,fj->p", Ef, Kf, Ef[:, :, p_in]))


def synthetic_K(f, x, y, xs, ys, v):
    """Ideal response of a single point scatterer at (xs, ys): K_ij = exp(-i w (d_i + d_j)/v)."""
    d = np.hypot(x - xs, y - ys)
    ph = np.exp(-2j * np.pi * f[:, None] * d[None] / v)          # [f, i]
    return ph[:, :, None] * ph[:, None, :]


# ---------------------------------------------------------------- diagnostics
def pair_snr(C, lags, dist, v_lo=1.0, v_hi=3.5, gap_s=10.0, noise_len_s=30.0):
    """Symmetric-lag SNR: max |C| inside d/v_hi..d/v_lo over the std of the noise window that
    starts gap_s after the signal window and lasts noise_len_s (kept close to the signal so the
    window-taper decay of the correlation does not shrink the noise estimate).
    Independent white noise gives SNR of about 3-4."""
    out = {}
    mid = len(lags) // 2                                  # index of zero lag
    m = np.arange(1, mid)                                 # positive lag samples
    t = lags[mid + m]
    for i in range(C.shape[0]):
        for j in range(i + 1, C.shape[0]):
            c = C[i, j]
            sym = 0.5 * (c[mid + m] + c[mid - m])         # average of +tau and -tau
            sig = (t >= dist[i, j] / v_hi) & (t <= dist[i, j] / v_lo)
            t_n0 = dist[i, j] / v_lo + gap_s
            noi = (t >= t_n0) & (t <= t_n0 + noise_len_s)
            out[(i, j)] = np.abs(sym[sig]).max() / (sym[noi].std() + 1e-30) if sig.any() else np.nan
    return out


# ---------------------------------------------------------------- driver
def run(st, coords=None, fmin=0.1, fmax=1.0, v=2.0, win_s=120.0, overlap=0.5,
        v_mute=1.0, t_pad=3.0, t_max=60.0, dx_km=0.5, pad_km=10.0, outdir=".", plot=True):
    names, fs, data, lat, lon = stream_to_array(st, coords)
    return run_array(names, fs, data, lat, lon, fmin, fmax, v, win_s, overlap,
                     v_mute, t_pad, t_max, dx_km, pad_km, outdir, plot)


def run_array(names, fs, data, lat, lon, fmin=0.1, fmax=1.0, v=2.0, win_s=120.0, overlap=0.5,
              v_mute=1.0, t_pad=3.0, t_max=60.0, dx_km=0.5, pad_km=10.0, outdir=".", plot=True):
    import os
    x, y, origin = to_local_km(lat, lon)
    dist = distance_matrix(x, y)

    print(f"[1] {len(names)} stations {names}, fs={fs} Hz, {data.shape[1]/fs/60:.1f} min")
    d = preprocess(data, fs, fmin, fmax)
    U, f, nfft = window_spectra(d, fs, win_s, overlap, fmin, fmax)
    print(f"[3] {U.shape[0]} windows of {win_s:.0f} s, df={f[1]-f[0]:.4f} Hz")
    K = response_matrix(U)
    C, lags = correlations(K, nfft, fs)
    snr = pair_snr(C, lags, dist)
    print("[4] correlation SNR per pair (symmetric lag, signal d/3.5..d/1.0 km/s; below ~6 = no signal):")
    for (i, j), s in sorted(snr.items(), key=lambda kv: dist[kv[0]]):
        print(f"      {names[i]}-{names[j]}  {dist[i, j]:5.1f} km  SNR {s:5.1f}")

    Cg = gate(C, lags, dist, v_mute, t_pad, t_max)
    Kg = correlations_to_K(Cg, nfft)
    fmask = (f >= fmin) & (f <= fmax)

    gx, gy, px, py = make_grid(x, y, pad_km, dx_km)
    E = steering(f, x, y, px, py, v)
    img = confocal_image(Kg, E, fmask).reshape(len(gy), len(gx))
    pk = np.argmax(img)
    row = focused_row(Kg, E, fmask, pk).reshape(img.shape)
    # PSF: what a perfect point scatterer at the image peak would look like with these 5 stations
    Ks = synthetic_K(f, x, y, px[pk], py[pk], v)
    psf = confocal_image(Ks, E, fmask).reshape(img.shape)
    print(f"[6] confocal image peak at x={px[pk]:.1f} km, y={py[pk]:.1f} km (v={v} km/s)")

    np.savez(os.path.join(outdir, "mi_results.npz"), names=names, x=x, y=y, origin=origin,
             f=f, K=K, Kg=Kg, lags=lags, C=C, Cg=Cg, gx=gx, gy=gy, img=img, row=row, psf=psf,
             params=dict(fmin=fmin, fmax=fmax, v=v, win_s=win_s, overlap=overlap,
                         v_mute=v_mute, t_pad=t_pad, t_max=t_max))
    if plot:
        _plots(names, x, y, dist, C, Cg, lags, gx, gy, img, row, psf, px[pk], py[pk], t_max, outdir)
    return dict(names=names, x=x, y=y, f=f, K=K, Kg=Kg, C=C, Cg=Cg, lags=lags,
                gx=gx, gy=gy, img=img, row=row, psf=psf, snr=snr)


def _plots(names, x, y, dist, C, Cg, lags, gx, gy, img, row, psf, xp, yp, t_max, outdir):
    import os
    import matplotlib.pyplot as plt
    n = len(names)
    pairs = sorted([(i, j) for i in range(n) for j in range(i + 1, n)], key=lambda p: dist[p])
    sel = np.abs(lags) <= t_max

    fig, axs = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    for ax, CC, title in [(axs[0], C, "Raw correlations"), (axs[1], Cg, "Direct wave muted")]:
        for (i, j) in pairs:
            c = CC[i, j, sel]
            c = c / (np.abs(C[i, j, sel]).max() + 1e-30)
            ax.plot(lags[sel], dist[i, j] + 1.5 * c, "k", lw=0.6)
            ax.text(t_max * 1.01, dist[i, j], f"{names[i]}-{names[j]}", fontsize=7, va="center")
        ax.set_xlabel("Lag (s)")
        ax.set_title(title)
    axs[0].set_ylabel("Inter-station distance (km)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "mi_correlations.png"), dpi=150)
    plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
    ext = [gx[0], gx[-1], gy[0], gy[-1]]
    for ax, A, title in [(axs[0], img, "Confocal image |R(r,r)|"),
                         (axs[1], row, "Focused matrix row |R(r_peak, r_out)|"),
                         (axs[2], psf, "Array PSF: synthetic scatterer at peak")]:
        im = ax.imshow(A / A.max(), origin="lower", extent=ext, cmap="viridis", vmin=0, vmax=1)
        ax.plot(x, y, "w^", ms=8, mec="k")
        for s, xi, yi in zip(names, x, y):
            ax.text(xi, yi + 0.8, s, color="w", fontsize=7, ha="center")
        ax.plot(xp, yp, "r+", ms=10)
        ax.set_xlabel("East (km)")
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, shrink=0.8, label="Normalised amplitude")
    axs[0].set_ylabel("North (km)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "mi_images.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    from obspy import read
    run(read("NZ.*.EHZ.sac"))
