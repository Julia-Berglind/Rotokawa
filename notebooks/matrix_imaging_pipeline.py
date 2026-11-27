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
  6. focus         project K onto a grid of virtual sources/receivers using travel times T_i(r, f)
                   through the Rotokawa velocity model (data/velocity_model.npz, notebook 01):
                     R(r_in, r_out) = sum_f sum_ij exp(+i w T_i(r_out, f)) K_ij exp(+i w T_j(r_in, f))
                   T_i is a 2D eikonal (fast marching) travel time from station i through a map of local
                   Rayleigh phase velocity c(x, y, f), estimated from each Vs column of the model as
                   c = RAYLEIGH_VS_RATIO * mean Vs over the top lambda/LAMBDA_FRAC (a rule of thumb, not a
                   dispersion calculation). Pass v=<km/s> with vmodel=None for the old constant velocity.
                   confocal image I(r) = |R(r, r)|, plus the focused matrix around one point
                   and the array point-spread function (PSF) for a synthetic point scatterer.

Coordinates are km east and north of the drill site in data/site.json, as in the other notebooks, so the
images line up with the velocity model and stations.csv.

Usage (in the rotokawa environment): see geonet_matrix_imaging.py, or
    import matrix_imaging_pipeline as mip
    out = mip.run(st)                  # writes figures + mi_results.npz to data/geonet_matrix_imaging/

Parameters that change the result (all defaults are my picks, not established values):
    fmin, fmax   correlation/imaging band (Hz)            default 0.1-1.0
    vmodel       velocity model .npz                       default data/velocity_model.npz
    v            constant velocity (km/s), only used when vmodel=None
    win_s        correlation window length (s)            default 120
    v_mute       direct-wave mute velocity (km/s)          default 1.0  (mutes |tau| < d/v_mute + t_pad)
"""
import json
import os
from pathlib import Path

import numpy as np
from scipy.signal import butter, sosfiltfilt, detrend
from scipy.ndimage import uniform_filter1d
from scipy.interpolate import RegularGridInterpolator

DATA = Path(__file__).resolve().parent.parent / "data"      # repo data/ folder (this file lives in notebooks/)
OUTDIR = DATA / "geonet_matrix_imaging"
VMODEL = DATA / "velocity_model.npz"
SITE = DATA / "site.json"
KY = 111.32                       # km per degree latitude (as in notebooks 01 and 02)
RAYLEIGH_VS_RATIO = 0.92          # Rayleigh phase velocity / Vs of a half-space (Poisson solid)
LAMBDA_FRAC = 3.0                 # Rayleigh waves sample roughly the top wavelength / 3
N_REF_FREQ = 8                    # frequencies at which travel times are computed (interpolated between)


# ---------------------------------------------------------------- 1. load
def stream_to_array(st, coords=None):
    """coords: optional {station: (lat, lon)}; otherwise read from tr.stats.sac.stla/stlo."""
    st = st.copy()
    st.sort(keys=["station"])
    # one vertical channel per station: prefer HHZ over EHZ if both came back
    keep = {}
    for tr in st:
        s = tr.stats.station
        if s not in keep or (tr.stats.channel.startswith("HH") and not keep[s].stats.channel.startswith("HH")):
            keep[s] = tr
    st.traces = [keep[s] for s in sorted(keep)]
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


def drill_site(site=SITE):
    s = json.load(open(site))["drill_site"]
    return s["lat"], s["lon"]


def to_local_km(lat, lon, lat0=None, lon0=None):
    """km east and north of (lat0, lon0); defaults to the drill site in data/site.json."""
    if lat0 is None or lon0 is None:
        lat0, lon0 = drill_site()
    x = (np.asarray(lon) - lon0) * np.cos(np.radians(lat0)) * KY
    y = (np.asarray(lat) - lat0) * KY
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
def make_grid(x, y, pad_km=10.0, dx_km=0.5, bounds=None):
    """Imaging grid around the stations; bounds=(xmin, xmax, ymin, ymax) clips it (e.g. to the velocity model)."""
    x0, x1, y0, y1 = x.min() - pad_km, x.max() + pad_km, y.min() - pad_km, y.max() + pad_km
    if bounds is not None:
        x0, x1, y0, y1 = max(x0, bounds[0]), min(x1, bounds[1]), max(y0, bounds[2]), min(y1, bounds[3])
    gx = np.arange(x0, x1 + dx_km / 2, dx_km)
    gy = np.arange(y0, y1 + dx_km / 2, dx_km)
    GX, GY = np.meshgrid(gx, gy)
    return gx, gy, GX.ravel(), GY.ravel()


# ---------------------------------------------------------------- 6a. velocity model -> travel times
def load_velocity_model(path=VMODEL):
    """data/velocity_model.npz from notebook 01: x, y, z (km), vs[x, y, z] (km/s)."""
    m = np.load(path)
    return dict(x=m["x"], y=m["y"], z=m["z"], vs=m["vs"], vp=m["vp"])


def rayleigh_phase_velocity(vm, freqs, ratio=RAYLEIGH_VS_RATIO, lam_frac=LAMBDA_FRAC, n_iter=10):
    """c[k, x, y]: approximate Rayleigh phase velocity at freqs[k] from each Vs column:
    c = ratio * mean(Vs over 0..lambda/lam_frac), lambda = c/f, solved by fixed-point iteration.
    Long wavelengths that reach the model base average the whole column."""
    z, vs = vm["z"], vm["vs"]
    dz = np.gradient(z)
    cum = np.cumsum(vs * dz, axis=-1)                        # integral of Vs from the surface
    cumz = np.cumsum(dz)
    out = np.empty((len(freqs),) + vs.shape[:2])
    for k, fk in enumerate(freqs):
        c = ratio * vs[..., 0].copy()
        for _ in range(n_iter):
            h = np.clip(c / fk / lam_frac, cumz[0], cumz[-1])  # sampling depth for each column
            iz = np.searchsorted(cumz, h).clip(0, len(z) - 1)
            mean_vs = np.take_along_axis(cum, iz[..., None], -1)[..., 0] / cumz[iz]
            c = ratio * mean_vs
        out[k] = c
    return out


def eikonal_2d(gx, gy, c, sx, sy, r0=None):
    """First-arrival time T[x, y] from (sx, sy) through the phase-velocity map c[x, y] (fast marching)."""
    import skfmm
    dx = [gx[1] - gx[0], gy[1] - gy[0]]
    r0 = 1.5 * max(dx) if r0 is None else r0
    X, Y = np.meshgrid(gx, gy, indexing="ij")
    d = np.hypot(X - sx, Y - sy)
    v0 = RegularGridInterpolator((gx, gy), c)([[sx, sy]])[0]
    t = skfmm.travel_time(d - r0, c, dx=dx, order=2)
    return np.where(d < r0, d / v0, np.asarray(t) + r0 / v0)


def travel_times(vm, f, fmin, fmax, x, y, px, py, n_ref=N_REF_FREQ):
    """T[f, i, p]: Rayleigh travel time from station i to pixel p at every frequency in f.
    Computed at n_ref log-spaced frequencies in [fmin, fmax] and interpolated linearly in between
    (constant outside the band). Also returns the reference frequencies and phase-velocity maps."""
    fr = np.geomspace(fmin, fmax, n_ref)
    cmap = rayleigh_phase_velocity(vm, fr)
    interp_pts = np.column_stack([px, py])
    Tref = np.empty((n_ref, len(x), len(px)))
    for k in range(n_ref):
        for i in range(len(x)):
            T = eikonal_2d(vm["x"], vm["y"], cmap[k], x[i], y[i])
            Tref[k, i] = RegularGridInterpolator((vm["x"], vm["y"]), T)(interp_pts)
    fc = np.clip(f, fr[0], fr[-1])
    j = np.clip(np.searchsorted(fr, fc) - 1, 0, n_ref - 2)
    w = ((fc - fr[j]) / (fr[j + 1] - fr[j]))[:, None, None]
    return (1 - w) * Tref[j] + w * Tref[j + 1], fr, cmap


def steering(f, T):
    """E[f, i, p] = exp(+i 2 pi f T_ip(f)): undoes the travel-time phase from station i to pixel p."""
    return np.exp(2j * np.pi * f[:, None, None] * T)


def constant_velocity_times(f, x, y, px, py, v):
    d = np.hypot(x[:, None] - px[None, :], y[:, None] - py[None, :])  # [i, p]
    return np.broadcast_to(d / v, (len(f),) + d.shape)


def confocal_image(K, E, fmask):
    """I(p) = | sum_f sum_ij E_ip K_ij E_jp |  (virtual source = virtual receiver = p)."""
    Kf, Ef = K[fmask], E[fmask]
    return np.abs(np.einsum("fip,fij,fjp->p", Ef, Kf, Ef))


def focused_row(K, E, fmask, p_in):
    """|R(r_in = p_in, r_out = every pixel)|: where energy focused at p_in comes back to."""
    Kf, Ef = K[fmask], E[fmask]
    return np.abs(np.einsum("fip,fij,fj->p", Ef, Kf, Ef[:, :, p_in]))


def synthetic_K(f, Ts):
    """Ideal response of a single point scatterer with station travel times Ts[f, i]:
    K_ij = exp(-i w (T_i + T_j))."""
    ph = np.exp(-2j * np.pi * f[:, None] * Ts)                   # [f, i]
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
def run(st, coords=None, fmin=0.1, fmax=1.0, vmodel=VMODEL, v=2.0, win_s=120.0, overlap=0.5,
        v_mute=1.0, t_pad=3.0, t_max=60.0, dx_km=0.5, pad_km=10.0, outdir=OUTDIR, plot=True):
    names, fs, data, lat, lon = stream_to_array(st, coords)
    return run_array(names, fs, data, lat, lon, fmin, fmax, vmodel, v, win_s, overlap,
                     v_mute, t_pad, t_max, dx_km, pad_km, outdir, plot)


def run_array(names, fs, data, lat, lon, fmin=0.1, fmax=1.0, vmodel=VMODEL, v=2.0, win_s=120.0, overlap=0.5,
              v_mute=1.0, t_pad=3.0, t_max=60.0, dx_km=0.5, pad_km=10.0, outdir=OUTDIR, plot=True):
    """vmodel: path to the velocity model .npz (or a dict from load_velocity_model); None = constant v."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    x, y, origin = to_local_km(lat, lon)
    dist = distance_matrix(x, y)

    print(f"[1] {len(names)} stations {list(names)}, fs={fs} Hz, {data.shape[1]/fs/60:.1f} min")
    print(f"    km from drill site: " + ", ".join(f"{n} ({xi:+.1f}, {yi:+.1f})" for n, xi, yi in zip(names, x, y)))
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
    fb = f[fmask]

    if vmodel is not None:
        vm = load_velocity_model(vmodel) if not isinstance(vmodel, dict) else vmodel
        m = 2 * dx_km                                            # stay inside the model
        bounds = (vm["x"][0] + m, vm["x"][-1] - m, vm["y"][0] + m, vm["y"][-1] - m)
        if not (np.all((x > bounds[0]) & (x < bounds[1])) and np.all((y > bounds[2]) & (y < bounds[3]))):
            raise ValueError("a station lies outside the velocity model grid")
        gx, gy, px, py = make_grid(x, y, pad_km, dx_km, bounds)
        T, fref, cref = travel_times(vm, fb, fmin, fmax, x, y, px, py)
        cmean = cref.mean(axis=(1, 2))
        print("[6] Rayleigh phase velocity from the velocity model (map mean): "
              + ", ".join(f"{fr:.2f} Hz {c:.2f} km/s" for fr, c in zip(fref, cmean)))
        vel_info = dict(vmodel=str(vmodel) if not isinstance(vmodel, dict) else "dict",
                        f_ref=fref, c_ref_mean=cmean, rayleigh_vs_ratio=RAYLEIGH_VS_RATIO, lambda_frac=LAMBDA_FRAC)
    else:
        gx, gy, px, py = make_grid(x, y, pad_km, dx_km)
        T = constant_velocity_times(fb, x, y, px, py, v)
        fref = cref = None
        vel_info = dict(vmodel=None, v=v)
        print(f"[6] constant velocity {v} km/s")

    E = steering(fb, T)
    Kb, Ksel = Kg[fmask], np.ones(len(fb), bool)              # E and K now share the band-limited f axis
    img = confocal_image(Kb, E, Ksel).reshape(len(gy), len(gx))
    pk = np.argmax(img)
    row = focused_row(Kb, E, Ksel, pk).reshape(img.shape)
    # PSF: what a perfect point scatterer at the image peak would look like with these stations
    Ks = synthetic_K(fb, T[:, :, pk])
    psf = confocal_image(Ks, E, Ksel).reshape(img.shape)
    print(f"[6] confocal image peak at x={px[pk]:.1f} km, y={py[pk]:.1f} km east/north of the drill site")

    np.savez(outdir / "mi_results.npz", names=names, x=x, y=y, origin=origin,
             f=f, K=K, Kg=Kg, lags=lags, C=C, Cg=Cg, gx=gx, gy=gy, img=img, row=row, psf=psf,
             f_ref=fref if fref is not None else np.zeros(0),
             c_ref=cref if cref is not None else np.zeros(0),
             params=dict(fmin=fmin, fmax=fmax, win_s=win_s, overlap=overlap,
                         v_mute=v_mute, t_pad=t_pad, t_max=t_max, dx_km=dx_km, pad_km=pad_km, **vel_info))
    if plot:
        _plots(names, x, y, dist, C, Cg, lags, gx, gy, img, row, psf, px[pk], py[pk], t_max, outdir)
        if vmodel is not None:
            _plot_velocity(vm, fref, cref, names, x, y, outdir)
    print(f"[7] results in {outdir}")
    return dict(names=names, x=x, y=y, f=f, K=K, Kg=Kg, C=C, Cg=Cg, lags=lags,
                gx=gx, gy=gy, img=img, row=row, psf=psf, snr=snr, f_ref=fref, c_ref=cref)


def _plot_velocity(vm, fref, cref, names, x, y, outdir):
    import matplotlib.pyplot as plt
    ks = sorted({0, len(fref) // 2, len(fref) - 1})
    fig, axs = plt.subplots(1, len(ks), figsize=(5 * len(ks), 4.6), squeeze=False)
    for ax, k in zip(axs[0], ks):
        pc = ax.pcolormesh(vm["x"], vm["y"], cref[k].T, cmap="viridis", shading="auto")
        ax.plot(x, y, "w^", ms=8, mec="k")
        for s, xi, yi in zip(names, x, y):
            ax.text(xi, yi + 0.8, s, color="w", fontsize=7, ha="center")
        ax.set(aspect="equal", xlabel="East of drill site (km)", title=f"Rayleigh phase velocity, {fref[k]:.2f} Hz")
        fig.colorbar(pc, ax=ax, shrink=0.8, label="km/s")
    axs[0, 0].set_ylabel("North of drill site (km)")
    fig.tight_layout()
    fig.savefig(Path(outdir) / "mi_phase_velocity.png", dpi=150)
    plt.close(fig)


def _plots(names, x, y, dist, C, Cg, lags, gx, gy, img, row, psf, xp, yp, t_max, outdir):
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
        ax.set_xlabel("East of drill site (km)")
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, shrink=0.8, label="Normalised amplitude")
    axs[0].set_ylabel("North of drill site (km)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "mi_images.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    from obspy import read
    run(read(str(OUTDIR / "waveforms" / "NZ.*.sac")))
