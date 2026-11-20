"""Matrix imaging of the BDT: one forward model, one imaging step and one scoring method for 05, 05b, 06 and 07.

The Earth is the background model of 01 plus the low-velocity zone (LVZ) of one 04 scenario. The LVZ is turned into
Born scatterers, echoes are made between every pair of stations that record at the same time (reflection matrix R),
and R is focused back into a confocal image on a 3D grid. The image is scored against the true BDT top.

    from mi_common import MatrixImaging
    mi = MatrixImaging(phase="p")           # loads data/ once
    mi.run(results="../data/results_05_mi.csv", image_dir="../data/mi_images")

Hooks used by the later notebooks:
    drop=...        velocity drop inside the LVZ (05b)
    floor=...       fixed detection floor instead of the image's own median (05b)
    noise=...       correlation noise added to R, shrinking with recording overlap (06)
    tt_scat=...     echo travel times through a different "true" model (07)
    ref="1d"        focus with the 1D tables instead of the 3D ones (07)
"""

# --- Imports ---
import json
import os
import time

import numpy as np
import pandas as pd
from scipy.interpolate import RectBivariateSpline

DATA = "../data"

# --- LVZ: keep in step with 04 ---
DROP_VP = 0.05      # fractional Vp drop inside the BDT zone
VS_RATIO = 2.0      # Vs drop / Vp drop (docs/02_vs_contrast_notes.md)
EDGE = 0.3          # km over which the drop ramps in at the top and base


class MatrixImaging:
    """Inputs, settings and functions for matrix imaging with P or S waves."""

    def __init__(self, phase="p", data=DATA, seed=4, **settings):
        # --- Settings (any can be overridden by keyword) ---
        self.PHASE = phase
        self.DATA = data
        self.DF = 0.1 if phase == "p" else 0.05   # Hz. small enough that two-way times (P < ~8 s, S < ~14 s) do not wrap
        self.F_BAND = (1.0, 6.0)                  # Hz
        self.DROP = DROP_VP if phase == "p" else VS_RATIO * DROP_VP   # drop in the velocity this run images with
        self.EDGE = EDGE
        self.BOX, self.BOX_TAPER = 9.0, 2.0       # scatterers within +-BOX km of the drill site, tapered over BOX_TAPER
        self.SC_STEP = (1, 1, 1)                  # scatterer grid = travel-time grid every (nx, ny, nz) nodes
        self.IM_R, self.IM_Z = 6.0, (2.0, 9.0)    # image within +-IM_R km of the drill site, IM_Z km deep
        self.IM_STEP = (2, 2, 1)                  # image grid = travel-time grid every (nx, ny, nz) nodes
        self.TARGET_R = 3.0                       # km: columns this close to the scenario centre count as the target
        self.SNR_MIN = 3.0                        # a column is imaged if its peak is SNR_MIN x the background
        self.CREST_FRAC = 0.5                     # BDT high picked only from columns this bright relative to the brightest
        self.TOP_TOL = 0.3                        # km: an imaged column covers the top if within this of the truth
        self.P2_SEARCH_R = 4.0                    # km: phase 2 centre kept this close to the planned centre (as in 02)
        self.DAYS_REF = 30                        # noise levels refer to this many days of shared recording (06)
        for k, v in settings.items():
            if not hasattr(self, k):
                raise ValueError(f"unknown setting {k}")
            setattr(self, k, v)
        self.F = np.arange(self.F_BAND[0], self.F_BAND[1] + 1e-9, self.DF)
        self.WF = np.hanning(len(self.F) + 2)[1:-1]   # source spectrum: Hann taper over the band
        self.rng = np.random.default_rng(seed)
        self._load()

    # ------------------------------------------------------------------------------------------------------------
    # Inputs and grids
    # ------------------------------------------------------------------------------------------------------------
    def _load(self):
        d = self.DATA
        vm = np.load(f"{d}/velocity_model.npz")
        self.vx, self.vy, self.vz = vm["x"], vm["y"], vm["z"]
        self.v = vm["vp"] if self.PHASE == "p" else vm["vs"]
        self.stations = pd.read_csv(f"{d}/stations.csv", index_col="station")
        self.D = json.load(open(f"{d}/designs.json"))
        self.sc_npz = np.load(f"{d}/scenarios.npz")
        self.SC = json.load(open(f"{d}/scenarios.json"))
        tg = np.load(f"{d}/tt/grid.npz")
        self.tx, self.ty, self.tz = tx, ty, tz = tg["x"], tg["y"], tg["z"]
        t1 = np.load(f"{d}/tt/tt1d.npz")
        self.tt1 = RectBivariateSpline(t1["r"], t1["z"], t1[f"t_{self.PHASE}"])

        # travel-time grid -> velocity-model grid indices (the tt grid is a crop of the model grid)
        self.mx, self.my, self.mz = (np.searchsorted(np.round(a, 6), np.round(b, 6))
                                     for a, b in ((self.vx, tx), (self.vy, ty), (self.vz, tz)))
        assert np.allclose(self.vx[self.mx], tx) and np.allclose(self.vy[self.my], ty) \
            and np.allclose(self.vz[self.mz], tz)

        # scatterer candidates: every node at depths any scenario's LVZ reaches, tapered towards the box edge
        f = self.sc_npz.files
        self.SC_Z = (np.nanmin([self.sc_npz[k] for k in f if k.startswith("top_")]) - self.EDGE,
                     np.nanmax([self.sc_npz[k] for k in f if k.startswith("base_")]) + self.EDGE)
        si, sj, sk = self.nodes(self.BOX, self.SC_Z, self.SC_STEP)
        self.SI, self.SJ, self.SK = (a.ravel() for a in np.meshgrid(si, sj, sk, indexing="ij"))
        edge = np.clip((self.BOX - np.maximum(np.abs(tx[self.SI]), np.abs(ty[self.SJ]))) / self.BOX_TAPER, 0, 1)
        self.SC_TAPER = 0.5 - 0.5 * np.cos(np.pi * edge)
        st = self.SC_STEP
        self.CELL = (tx[1] - tx[0]) * st[0] * (ty[1] - ty[0]) * st[1] * (tz[1] - tz[0]) * st[2]

        # image grid
        self.ii, self.ij, self.ik = self.nodes(self.IM_R, self.IM_Z, self.IM_STEP)
        self.IX, self.IY, self.IZ = tx[self.ii], ty[self.ij], tz[self.ik]
        self.II, self.IJ, self.IK = (a.ravel() for a in np.meshgrid(self.ii, self.ij, self.ik, indexing="ij"))
        self.XX, self.YY = np.meshgrid(self.IX, self.IY, indexing="ij")

        # deployment
        self.DEPLOY_END = max(p["end"] for ph in self.D["designs"].values() for p in ph)
        self.PERMANENT = self.D["trillium"] + self.D.get("existing", [])
        self.P2 = self.D["p2"]
        self.AZR = np.radians(self.P2["az"])
        self.BANK = {(i, j): s for s, i, j in zip(self.D["bank"]["ids"], self.D["bank"]["i"], self.D["bank"]["j"])}
        self._tt, self._tt1 = {}, {}

    def __repr__(self):
        return (f"MatrixImaging({self.PHASE.upper()}): image grid {len(self.ii)} x {len(self.ij)} x {len(self.ik)}, "
                f"{len(self.SI)} scatterer candidates ({self.SC_Z[0]:.1f}-{self.SC_Z[1]:.1f} km deep), "
                f"{len(self.F)} frequencies")

    def nodes(self, xr, zr, step):
        """Travel-time-grid node indices within +-xr km horizontally and zr depth, every step nodes."""
        i = np.where(np.abs(self.tx) <= xr + 1e-9)[0][::step[0]]
        j = np.where(np.abs(self.ty) <= xr + 1e-9)[0][::step[1]]
        k = np.where((self.tz >= zr[0] - 1e-9) & (self.tz <= zr[1] + 1e-9))[0][::step[2]]
        return i, j, k

    def scenarios(self, shapes=None):
        """Scenario names, optionally only those with a shape in shapes."""
        return [n for n, m in self.SC["scenarios"].items() if shapes is None or m["shape"] in shapes]

    # ------------------------------------------------------------------------------------------------------------
    # Earth: the LVZ of one scenario as Born scatterers
    # ------------------------------------------------------------------------------------------------------------
    def perturbation(self, name, drop=None):
        """Scatterer nodes (indices into SI/SJ/SK) and strengths sigma = d(1/v^2) x cell volume x edge taper."""
        drop = self.DROP if drop is None else drop
        mx, my = self.mx[self.SI], self.my[self.SJ]
        top, base = self.sc_npz[f"top_{name}"][mx, my], self.sc_npz[f"base_{name}"][mx, my]
        zz = self.tz[self.SK]
        w = np.clip(np.minimum(zz - np.nan_to_num(top, nan=99), np.nan_to_num(base, nan=99) - zz) / self.EDGE + 0.5,
                    0, 1)
        v0 = self.v[mx, my, self.mz[self.SK]]
        sig = (1 / (v0 * (1 - drop * w)) ** 2 - 1 / v0 ** 2) * self.CELL * self.SC_TAPER
        keep = np.abs(sig) > 1e-6 * np.abs(sig).max()
        return np.where(keep)[0], sig[keep].astype(np.float32)

    # ------------------------------------------------------------------------------------------------------------
    # Travel times and spreading
    # ------------------------------------------------------------------------------------------------------------
    def station_tt(self, s):
        """(3D travel times to all scatterer candidates, to all image nodes) for station s, from the 03 tables."""
        if s not in self._tt:
            t = np.load(f"{self.DATA}/tt/{s:03d}_{self.PHASE}.npy", mmap_mode="r")
            self._tt[s] = (np.asarray(t[self.SI, self.SJ, self.SK]), np.asarray(t[self.II, self.IJ, self.IK]))
        return self._tt[s]

    def tt_image(self, s, ref="3d"):
        """Focusing travel times from station s to the image nodes: '3d' (03 tables) or '1d' (03 1D table)."""
        if ref == "3d":
            return self.station_tt(s)[1]
        if ref != "1d":
            raise ValueError(f"unknown focusing model {ref}")
        if s not in self._tt1:
            r = np.hypot(self.tx[self.II] - self.stations.x[s], self.ty[self.IJ] - self.stations.y[s])
            self._tt1[s] = self.tt1.ev(r, self.tz[self.IK]).astype(np.float32)
        return self._tt1[s]

    def spreading(self, s, x, y, z):
        r = np.sqrt((x - self.stations.x[s]) ** 2 + (y - self.stations.y[s]) ** 2 + z ** 2)
        return (1 / (4 * np.pi * np.maximum(r, 0.05))).astype(np.float32)

    # ------------------------------------------------------------------------------------------------------------
    # Forward model, noise and confocal image
    # ------------------------------------------------------------------------------------------------------------
    def forward(self, ids, k_sc, sig, tt_scat=None):
        """Reflection matrices R(f) between stations ids for scatterers k_sc. (nF, N, N)

        tt_scat(s) -> travel times from s to all scatterer candidates; default the 3D tables (the model is right).
        """
        tt_scat = tt_scat or (lambda s: self.station_tt(s)[0])
        T = np.stack([tt_scat(s)[k_sc] for s in ids]).astype(np.float32)
        A = np.stack([self.spreading(s, self.tx[self.SI[k_sc]], self.ty[self.SJ[k_sc]], self.tz[self.SK[k_sc]])
                      for s in ids])
        R = np.empty((len(self.F), len(ids), len(ids)), np.complex64)
        for n, f in enumerate(self.F):
            G = A * np.exp(-2j * np.pi * f * T)                      # G(k, s_i), (N, K)
            R[n] = (G * sig) @ G.T
        return R

    def add_noise(self, R, days, level):
        """R plus symmetric complex noise with std level x sqrt(DAYS_REF / days) for each station pair.

        days: (N, N) days each pair recorded together. Pairs that never did are masked out of the image anyway.
        """
        E = (self.rng.normal(size=R.shape) + 1j * self.rng.normal(size=R.shape)) / np.sqrt(2)
        E = (E + E.transpose(0, 2, 1)) / np.sqrt(2)
        std = level * np.sqrt(self.DAYS_REF / np.maximum(days, 1e-9)) * (days > 0)
        return (R + E * std).astype(np.complex64)

    def image(self, ids, R, pair_mask, ref="3d"):
        """Confocal image I(x) = |sum_f W(f) sum_ij G*(x, s_i) R_ij G*(s_j, x)| on the image grid."""
        T = np.stack([self.tt_image(s, ref) for s in ids]).T.astype(np.float32)   # (P, N)
        A = np.stack([self.spreading(s, self.tx[self.II], self.ty[self.IJ], self.tz[self.IK]) for s in ids]).T
        I = np.zeros(len(self.II), np.complex128)
        for n, f in enumerate(self.F):
            Gc = A * np.exp(2j * np.pi * f * T)
            I += self.WF[n] * np.einsum("pi,ij,pj->p", Gc, R[n] * pair_mask, Gc, optimize=True)
        return np.abs(I).reshape(len(self.ii), len(self.ij), len(self.ik)).astype(np.float32)

    # ------------------------------------------------------------------------------------------------------------
    # Who records with whom
    # ------------------------------------------------------------------------------------------------------------
    def windows(self, phases):
        """Recording windows per station: permanent stations all deployment, nodes their own phase."""
        w = {s: [(0, self.DEPLOY_END)] for s in self.PERMANENT}
        for p in phases:
            for s in p["stations"]:
                w.setdefault(s, []).append((p["start"], p["end"]))
        return w

    @staticmethod
    def overlap_days(ids, win):
        """Days stations i != j recorded at the same time. (N, N)"""
        M = np.zeros((len(ids), len(ids)))
        for a, s in enumerate(ids):
            for b in range(a + 1, len(ids)):
                M[a, b] = M[b, a] = sum(max(0, min(e1, e2) - max(s1, s2))
                                        for s1, e1 in win[s] for s2, e2 in win[ids[b]])
        return M

    def pair_mask(self, ids, win):
        """True where stations i != j record at the same time."""
        return self.overlap_days(ids, win) > 0

    # ------------------------------------------------------------------------------------------------------------
    # Scores
    # ------------------------------------------------------------------------------------------------------------
    def pick_top(self, img, min_peak=None):
        """Picked BDT top per image column: shallowest local maximum above half the column peak. NaN if not imaged."""
        peak = img.max(axis=2)
        imaged = peak >= (self.SNR_MIN * np.median(img) if min_peak is None else min_peak)
        a = img / np.maximum(peak[..., None], 1e-30)
        local = np.zeros_like(a, bool)
        local[..., 1:-1] = (a[..., 1:-1] >= a[..., :-2]) & (a[..., 1:-1] >= a[..., 2:]) & (a[..., 1:-1] >= 0.5)
        first = np.argmax(local, axis=2)
        return np.where(imaged & local.any(axis=2), self.IZ[first], np.nan)

    def crest(self, top, tol=0.1):
        """BDT high: centroid of the columns within tol km of the shallowest top (a plateau for flat-topped shapes)."""
        hi = np.isfinite(top) & (top <= np.nanmin(top) + tol)
        return np.array([self.XX[hi].mean(), self.YY[hi].mean()])

    def true_top(self, name):
        """True BDT top on the image map grid."""
        return self.sc_npz[f"top_{name}"][np.ix_(self.mx[self.ii], self.my[self.ij])]

    def score(self, img, name, floor=None):
        """Coverage of the BDT top, depth error, and how far the picked BDT high is from the true one.

        floor=None: a column is imaged if its peak is SNR_MIN x the image's own median (05).
        floor=value: ... SNR_MIN x this fixed floor instead, e.g. the background of a reference image (05b).
        """
        centre = np.array(self.SC["scenarios"][name]["centre"])
        true = self.true_top(name)
        est = self.pick_top(img, None if floor is None else self.SNR_MIN * floor)
        target = (np.hypot(self.XX - centre[0], self.YY - centre[1]) <= self.TARGET_R) & np.isfinite(true)
        err = np.abs(est - true)
        good = target & (err <= self.TOP_TOL)
        crest_true = self.crest(true)
        min_peak = self.CREST_FRAC * img.max() if floor is None else max(self.CREST_FRAC * img.max(),
                                                                          self.SNR_MIN * floor)
        strong = self.pick_top(img, min_peak)
        if np.isfinite(strong).any():
            crest_est = self.crest(strong)
            crest_err = float(np.hypot(*(crest_est - crest_true)))
        else:
            crest_est, crest_err = np.array([np.nan, np.nan]), np.nan
        hit = target & np.isfinite(est)
        return dict(coverage=good.sum() / max(target.sum(), 1),
                    top_err_km=float(np.nanmedian(err[hit])) if hit.any() else np.nan,
                    crest_err_km=crest_err, crest_x=crest_est[0], crest_y=crest_est[1])

    # ------------------------------------------------------------------------------------------------------------
    # Adaptive phase 2: centre the rolling patches on the BDT high picked from the phase 1 image
    # ------------------------------------------------------------------------------------------------------------
    def adaptive_phase2(self, crest_xy):
        """Phase 2 patches (15 x 8 bank sites, cut into 5 x 4 patches in the planned order) around crest_xy."""
        P2, AZR = self.P2, self.AZR
        c0 = np.array(P2["centre0"])
        d = np.asarray(crest_xy, float) - c0
        if np.hypot(*d) > self.P2_SEARCH_R:
            d *= self.P2_SEARCH_R / np.hypot(*d)
        along = d[0] * np.sin(AZR) + d[1] * np.cos(AZR)
        across = -d[0] * np.cos(AZR) + d[1] * np.sin(AZR)
        i0 = int(round(along / P2["spacing"])) - 7
        j0 = int(round(across / P2["spacing"] - 0.5)) - 3
        phases = []
        for q, (pi, pj) in enumerate(P2["order"]):
            ids = [self.BANK[(i0 + 5 * pi + a, j0 + 4 * pj + b)] for a in range(5) for b in range(4)
                   if (i0 + 5 * pi + a, j0 + 4 * pj + b) in self.BANK]
            t0 = P2["p1_days"] + q * P2["patch_days"]
            phases.append(dict(start=t0, end=t0 + P2["patch_days"], stations=ids))
        return phases, c0 + d

    # ------------------------------------------------------------------------------------------------------------
    # One scenario, every design
    # ------------------------------------------------------------------------------------------------------------
    def design_names(self, designs=None):
        """Designs to run: all in designs.json plus 'phased (adaptive)', or the subset asked for."""
        names = list(self.D["designs"]) + (["phased (adaptive)"] if "phased" in self.D["designs"] else [])
        return names if designs is None else [d for d in names if d in designs]

    def run_scenario(self, name, designs=None, drop=None, ref="3d", tt_scat=None, noise=None):
        """Image one scenario with each design. Returns (rows of scores, {design: image}).

        drop: LVZ velocity drop (default self.DROP). ref: focusing model, '3d' or '1d'.
        tt_scat: echo travel times through a different true model (see forward). noise: noise level (see add_noise).
        """
        todo = self.design_names(designs)
        k_sc, sig = self.perturbation(name, drop)
        static = [d for d in todo if d != "phased (adaptive)"]
        DS = self.D["designs"]
        need = sorted(set(self.PERMANENT) | {s for d in static for p in DS[d] for s in p["stations"]}
                      | (set(DS["phased"][0]["stations"]) if "phased (adaptive)" in todo else set()))
        R = self.forward(need, k_sc, sig, tt_scat)
        pos = {s: n for n, s in enumerate(need)}
        meta = self.SC["scenarios"][name]
        rows, imgs = [], {}

        def one(design, phases, Rsub, ids, extra=None):
            days = self.overlap_days(ids, self.windows(phases))
            if noise is not None:
                Rsub = self.add_noise(Rsub, days, noise)
            img = self.image(ids, Rsub, days > 0, ref)
            imgs[design] = img
            rows.append(dict(scenario=name, conductor=meta["conductor"], shape=meta["shape"], design=design,
                             ref=ref, n_stations=len(ids), n_pairs=int((days > 0).sum() // 2),
                             **self.score(img, name), **(extra or {})))
            return img

        def sub(ids):
            sel = [pos[s] for s in ids]
            return R[:, sel][:, :, sel]

        for d in static:
            ids = sorted(set(self.PERMANENT) | {s for p in DS[d] for s in p["stations"]})
            one(d, DS[d], sub(ids), ids)

        if "phased (adaptive)" in todo:
            p1 = DS["phased"][0]
            ids1 = sorted(set(self.PERMANENT) | set(p1["stations"]))
            pick = self.score(one("phased: phase 1 only", [p1], sub(ids1), ids1), name)
            rows.pop()                                   # phase 1 alone is not a design; keep only its pick
            imgs.pop("phased: phase 1 only")
            ok = np.isfinite(pick["crest_x"])
            xy = (pick["crest_x"], pick["crest_y"]) if ok else self.P2["centre0"]
            p2_phases, centre = self.adaptive_phase2(xy)
            phases = [p1] + p2_phases
            ids = sorted(set(self.PERMANENT) | {s for p in phases for s in p["stations"]})
            one("phased (adaptive)", phases, self.forward(ids, k_sc, sig, tt_scat), ids,
                extra=dict(p1_error_km=pick["crest_err_km"], p2_centre_x=centre[0], p2_centre_y=centre[1]))
        return rows, imgs

    # ------------------------------------------------------------------------------------------------------------
    # Every scenario x every design, saved as it goes
    # ------------------------------------------------------------------------------------------------------------
    def run(self, scenarios=None, designs=None, results=None, image_dir=None, extra=None, verbose=True, **kw):
        """Run run_scenario over scenarios and append the scores to results (CSV). Rerunning skips rows already there.

        extra: columns that identify this run in the results file, e.g. dict(nsr=0.5) or dict(variant="baseline").
        image_dir: save each scenario's images there (None = do not save). kw: passed to run_scenario.
        Returns the full results table.
        """
        extra = extra or {}
        names = self.scenarios() if scenarios is None else list(scenarios)
        wanted = self.design_names(designs)
        done = pd.read_csv(results) if results and os.path.exists(results) else pd.DataFrame()
        if image_dir:
            os.makedirs(image_dir, exist_ok=True)

        def have(name, d):
            if done.empty:
                return False
            m = (done.scenario == name) & (done.design == d)
            for k, v in extra.items():
                m &= (done[k] == v) if k in done else False
            return bool(m.any())

        for name in names:
            todo = [d for d in wanted if not have(name, d)]
            if not todo:
                continue
            t0 = time.time()
            rows, imgs = self.run_scenario(name, designs=todo, **kw)
            done = pd.concat([done, pd.DataFrame([{**r, **extra} for r in rows])], ignore_index=True)
            if results:
                done.to_csv(results, index=False)
            if image_dir:
                tag = "_".join(f"{v}" for v in extra.values())
                np.savez_compressed(f"{image_dir}/{name}{'_' + tag if tag else ''}.npz",
                                    x=self.IX, y=self.IY, z=self.IZ,
                                    **{k.replace(" ", "_"): v for k, v in imgs.items()})
            if verbose:
                print(f"{name}{' ' + str(extra) if extra else ''}: {', '.join(todo)} in {time.time() - t0:.0f} s",
                      flush=True)
        return done

    # ------------------------------------------------------------------------------------------------------------
    # Plot helper
    # ------------------------------------------------------------------------------------------------------------
    def section(self, img_file, name, half=6.0):
        """Images and true top along the node axis through the scenario centre.

        Returns (s, z, {design: section normalised to its max}, true top along s).
        """
        im = np.load(img_file)
        keys = [k for k in im.files if k not in ("x", "y", "z")]
        centre = self.SC["scenarios"][name]["centre"]
        s = np.linspace(-half, half, int(20 * half) + 1)
        px, py = centre[0] + s * np.sin(self.AZR), centre[1] + s * np.cos(self.AZR)
        jx = np.clip(np.round((px - im["x"][0]) / (im["x"][1] - im["x"][0])).astype(int), 0, len(im["x"]) - 1)
        jy = np.clip(np.round((py - im["y"][0]) / (im["y"][1] - im["y"][0])).astype(int), 0, len(im["y"]) - 1)
        true = self.sc_npz[f"top_{name}"][np.searchsorted(self.vx, im["x"][jx]), np.searchsorted(self.vy, im["y"][jy])]
        secs = {}
        for k in keys:
            sec = im[k][jx, jy, :]
            secs[k.replace("_", " ")] = sec / sec.max()
        return s, im["z"], secs, true
