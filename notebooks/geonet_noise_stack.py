"""
geonet_noise_stack.py
Quick look at how GeoNet noise correlations converge, and whether the direct surface wave matches the velocity
model. No imaging.

  1. Download N_HOURS one-hour blocks (in parallel), remove response, resample to 5 Hz, and reduce each hour to its
     response matrix K[f, i, j] (5x5 per frequency). Each hour is cached, so reruns take seconds.
  2. Convergence: stack 1, 2, 4, ... hours. Per pair, measure the correlation noise (std in a late-lag window where
     no direct wave arrives) and the direct-wave SNR, and fit noise ~ hours^-a. a = 0.5 is the sqrt(time) law that
     05_noise assumes; if it holds, extrapolate the noise to DAYS_REF days.
  3. Velocity check: in a few frequency bands, pick the envelope peak of each pair's full-stack correlation and
     compare with the group travel time predicted from the velocity model (via mip.rayleigh_phase_velocity).

In a notebook (from notebooks/):
    import geonet_noise_stack as gns
    res = gns.run()                       # or gns.run(n_hours=48, t0="2026-11-18T00:00:00")
Outputs: data/geonet_matrix_imaging/noise_stack/
           hours/K_<time>.npz          cached hourly response matrices
           convergence.csv, convergence.png, velocity_check.csv, velocity_check.png
"""
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

import matrix_imaging_pipeline as mip

STATIONS = ["ALRZ", "ARAZ", "THQ2", "THSZ", "WPRZ"]
CHANNELS = "HHZ,EHZ"
T0 = "2026-11-20T00:00:00"
N_HOURS = 24
FS = 5.0                                   # Hz after resampling (fmax is 1 Hz)
FMIN, FMAX = 0.1, 1.0
WIN_S, OVERLAP = 120.0, 0.5
PRE_FILT = (0.03, 0.05, 2.0, 2.4)          # Hz, response removal (applied before resampling)
NOISE_LAG = (70.0, 110.0)                  # s, late-lag window for the convergence curve (direct wave is < ~25 s);
                                           # the window taper shrinks it, but only the ratio between stacks is used
BANDS = [(0.1, 0.2), (0.2, 0.4), (0.4, 0.8)]
DAYS_REF = 30                              # as in 05_noise
WORKERS = 4
OUT = mip.OUTDIR / "noise_stack"


# ---------------------------------------------------------------- 1. hourly response matrices
def _hour_file(t):
    return OUT / "hours" / f"K_{t.strftime('%Y%m%dT%H')}.npz"


def hour_K(t, inv):
    """Download one hour from t, process, and save K. Returns the cache path, or None if a station is missing."""
    from obspy.clients.fdsn import Client
    path = _hour_file(t)
    if path.exists():
        return path
    st = Client("GEONET").get_waveforms("NZ", ",".join(STATIONS), "*", CHANNELS, t - 60, t + 3660)
    st.merge(fill_value=0)
    st.remove_response(inventory=inv, output="VEL", pre_filt=PRE_FILT)
    st.resample(FS).trim(t, t + 3600)
    for tr in st:
        co = inv.get_coordinates(tr.id, t)
        tr.stats.sac = {"stla": co["latitude"], "stlo": co["longitude"]}
    names, fs, data, lat, lon = mip.stream_to_array(st)
    if names != STATIONS or data.shape[1] < 0.9 * 3600 * FS:
        return None
    d = mip.preprocess(data, fs, FMIN, FMAX)
    U, f, nfft = mip.window_spectra(d, fs, WIN_S, OVERLAP, FMIN, FMAX)
    np.savez(path, K=mip.response_matrix(U), f=f, nfft=nfft, fs=fs, lat=lat, lon=lon, n_win=U.shape[0])
    return path


def download(t0=T0, n_hours=N_HOURS):
    from obspy import UTCDateTime
    from obspy.clients.fdsn import Client
    (OUT / "hours").mkdir(parents=True, exist_ok=True)
    t0 = UTCDateTime(t0)
    hours = [t0 + 3600 * k for k in range(n_hours)]
    todo = [t for t in hours if not _hour_file(t).exists()]
    if todo:
        inv = Client("GEONET").get_stations(network="NZ", station=",".join(STATIONS), channel=CHANNELS,
                                            level="response", starttime=hours[0], endtime=hours[-1] + 3600)
        print(f"downloading {len(todo)} hours with {WORKERS} workers ...")
        with ThreadPoolExecutor(WORKERS) as ex:
            jobs = {ex.submit(hour_K, t, inv): t for t in todo}
            for k, job in enumerate(as_completed(jobs), 1):
                try:
                    ok = job.result() is not None
                except Exception as e:                      # one bad hour should not stop the run
                    ok = False
                    print(f"  {jobs[job]}: {e}")
                print(f"  {k}/{len(todo)} {jobs[job].strftime('%m-%d %H:00')} {'ok' if ok else 'skipped'}",
                      flush=True)
    files = [_hour_file(t) for t in hours if _hour_file(t).exists()]
    print(f"{len(files)}/{n_hours} usable hours")
    return files


# ---------------------------------------------------------------- 2. convergence
def _sym(C, lags):
    mid = len(lags) // 2
    m = np.arange(1, mid)
    return lags[mid + m], 0.5 * (C[..., mid + m] + C[..., mid - m])


def convergence(Ks, nfft, fs, dist):
    """Noise std and SNR per pair for stacks of 1, 2, 4, ... hours (first n hours)."""
    n_list = sorted({2 ** k for k in range(int(np.log2(len(Ks))) + 1)} | {len(Ks)})
    rows = []
    for n in n_list:
        C, lags = mip.correlations(np.mean(Ks[:n], axis=0), nfft, fs)
        t, sym = _sym(C, lags)
        late = (t >= NOISE_LAG[0]) & (t <= NOISE_LAG[1])
        for i in range(len(dist)):
            for j in range(i + 1, len(dist)):
                sig = (t >= dist[i, j] / 3.5) & (t <= dist[i, j] / 1.0)
                near = (t >= dist[i, j] + 10) & (t <= dist[i, j] + 40)     # as mip.pair_snr
                rows.append(dict(hours=n, pair=f"{STATIONS[i]}-{STATIONS[j]}", dist_km=dist[i, j],
                                 noise=sym[i, j, late].std(),
                                 snr=np.abs(sym[i, j, sig]).max() / sym[i, j, near].std()))
    df = pd.DataFrame(rows)
    df["noise_rel"] = df.noise / df.groupby("pair").noise.transform("first")   # relative to 1 hour
    return df


def fit_exponent(df):
    """noise ~ hours^-a, fitted to the median over pairs."""
    med = df.groupby("hours").noise_rel.median()
    a = -np.polyfit(np.log(med.index), np.log(med.values), 1)[0]
    return a, med


# ---------------------------------------------------------------- 3. velocity check
def bandpass_K(K, f, lo, hi):
    return K * mip.band_taper(f, lo, hi, frac=0.2)[:, None, None]


def predicted_group_times(vm, x, y, bands, n_pts=60):
    """Group travel time per pair and band: straight path through the group-velocity map
    U = c / (1 - (f/c) dc/df), from the model's Rayleigh phase-velocity maps."""
    from scipy.interpolate import RegularGridInterpolator
    out = {}
    for lo, hi in bands:
        fc = np.sqrt(lo * hi)
        fr = fc * np.array([0.9, 1.0, 1.1])
        c = mip.rayleigh_phase_velocity(vm, fr)
        dcdf = (c[2] - c[0]) / (fr[2] - fr[0])
        U = c[1] / (1 - fc / c[1] * dcdf)
        interp = RegularGridInterpolator((vm["x"], vm["y"]), 1 / U)
        for i in range(len(x)):
            for j in range(i + 1, len(x)):
                s = np.linspace(0, 1, n_pts)
                px, py = x[i] + s * (x[j] - x[i]), y[i] + s * (y[j] - y[i])
                out[(i, j, lo, hi)] = interp(np.column_stack([px, py])).mean() * np.hypot(x[j] - x[i], y[j] - y[i])
    return out


def velocity_check(K, f, nfft, fs, x, y, dist, vm):
    from scipy.signal import hilbert
    pred = predicted_group_times(vm, x, y, BANDS)
    rows, traces = [], {}
    for lo, hi in BANDS:
        C, lags = mip.correlations(bandpass_K(K, f, lo, hi), nfft, fs)
        t, sym = _sym(C, lags)
        env = np.abs(hilbert(sym, axis=-1))
        for i in range(len(x)):
            for j in range(i + 1, len(x)):
                win = (t >= dist[i, j] / 4.0) & (t <= dist[i, j] / 0.7 + 3)
                k = np.argmax(np.where(win, env[i, j], 0))
                rows.append(dict(band=f"{lo}-{hi} Hz", pair=f"{STATIONS[i]}-{STATIONS[j]}", dist_km=dist[i, j],
                                 t_obs=t[k], t_pred=pred[(i, j, lo, hi)],
                                 u_obs=dist[i, j] / t[k], u_pred=dist[i, j] / pred[(i, j, lo, hi)]))
                traces[(lo, hi, i, j)] = (t, sym[i, j], env[i, j])
    return pd.DataFrame(rows), traces


# ---------------------------------------------------------------- plots
def plot_convergence(df, a, med, path):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for pair, g in df.groupby("pair"):
        axs[0].loglog(g.hours, g.noise_rel, color="0.75", lw=0.8)
        axs[1].semilogx(g.hours, g.snr, color="0.75", lw=0.8)
    h = med.index.values
    axs[0].loglog(h, med.values, "ko-", label=f"median over pairs: hours^-{a:.2f}")
    axs[0].loglog(h, h ** -0.5, "--", color="#d1495b", label="sqrt(time) law (05_noise)")
    axs[0].set(xlabel="Hours stacked", ylabel="Correlation noise / 1-hour noise", title="Noise convergence")
    axs[0].legend(fontsize=8)
    axs[1].semilogx(h, df.groupby("hours").snr.median().values, "ko-", label="median over pairs")
    axs[1].axhline(6, color="#d1495b", ls="--", lw=1, label="~detection")
    axs[1].set(xlabel="Hours stacked", ylabel="Direct-wave SNR", title="Direct-wave SNR")
    axs[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    return fig


def plot_velocity(vdf, traces, dist, path):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, len(BANDS), figsize=(5 * len(BANDS), 4.6), sharey=True)
    for ax, (lo, hi) in zip(np.atleast_1d(axs), BANDS):
        for (b0, b1, i, j), (t, s, e) in traces.items():
            if (b0, b1) != (lo, hi):
                continue
            sc = 1.5 / (np.abs(s[t < 60]).max() + 1e-30)
            ax.plot(t[t < 60], dist[i, j] + sc * s[t < 60], "k", lw=0.5)
        g = vdf[vdf.band == f"{lo}-{hi} Hz"]
        ax.plot(g.t_obs, g.dist_km, "o", color="#2a78d6", label="picked envelope peak")
        ax.plot(g.t_pred, g.dist_km, "x", color="#d1495b", ms=8, mew=2, label="predicted from velocity model")
        ax.set(xlabel="Lag (s, symmetric)", title=f"{lo}-{hi} Hz")
    np.atleast_1d(axs)[0].set_ylabel("Inter-station distance (km)")
    np.atleast_1d(axs)[0].legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    return fig


# ---------------------------------------------------------------- driver
def run(t0=T0, n_hours=N_HOURS, vmodel=mip.VMODEL, plot=True):
    files = download(t0, n_hours)
    if not files:
        raise RuntimeError("no usable hours")
    hs = [np.load(p) for p in files]
    Ks = np.array([h["K"] for h in hs])
    f, nfft, fs = hs[0]["f"], int(hs[0]["nfft"]), float(hs[0]["fs"])
    x, y, _ = mip.to_local_km(hs[0]["lat"], hs[0]["lon"])
    dist = mip.distance_matrix(x, y)

    cdf = convergence(Ks, nfft, fs, dist)
    a, med = fit_exponent(cdf)
    cdf.to_csv(OUT / "convergence.csv", index=False)
    n = len(Ks)
    days_ratio = DAYS_REF * 24 / n
    print(f"\nnoise ~ hours^-{a:.2f} (sqrt law: 0.50), median over pairs, {n} hours")
    print(f"extrapolated to {DAYS_REF} days: noise x{days_ratio ** -a:.3f} of the {n}-hour stack "
          f"(sqrt law would give x{days_ratio ** -0.5:.3f})")
    print("direct-wave SNR at the full stack:")
    print(cdf[cdf.hours == n].sort_values("dist_km")[["pair", "dist_km", "snr"]].round(1).to_string(index=False))

    vm = mip.load_velocity_model(vmodel)
    vdf, traces = velocity_check(Ks.mean(axis=0), f, nfft, fs, x, y, dist, vm)
    vdf.to_csv(OUT / "velocity_check.csv", index=False)
    print("\ngroup velocity, observed vs velocity model (median over pairs, km/s):")
    print(vdf.groupby("band")[["u_obs", "u_pred"]].median().round(2).to_string())
    print("(only trust pairs whose SNR above is > ~6)")

    if plot:
        plot_convergence(cdf, a, med, OUT / "convergence.png")
        plot_velocity(vdf, traces, dist, OUT / "velocity_check.png")
    print(f"\nresults in {OUT}")
    return dict(convergence=cdf, exponent=a, velocity=vdf, K=Ks.mean(axis=0), f=f, nfft=nfft, fs=fs,
                names=STATIONS, x=x, y=y)


if __name__ == "__main__":
    run()
