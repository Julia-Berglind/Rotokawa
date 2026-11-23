"""Matrix imaging of the BDT, shared by 05 (fully converged correlations) and 06 (with noise).

The Earth is the velocity model of 01 with a dome-shaped BDT zone in which Vp is DVP higher. The zone is turned into
Born scatterers, echoes are made between every pair of stations that record at the same time (the reflection
matrix R: what fully converged ambient-noise correlations give), and R is focused back into a confocal image with
the 3D travel times of 03.

Phased designs are run the way they would be in the field: phase 1 is imaged on its own, the BDT high is picked
from that image, and phase 2 is laid out around it with the same rules as 02 (patches on the station bank, sites in
water replaced by the nearest free dry site).

    from mi_common import MatrixImaging
    mi = MatrixImaging()                               # dome under the drill site
    out = mi.run_design("phased: interleaved", "dome")
    out["image"], out["phases"], out["score"]

06 adds correlation noise with run_design(..., noise=mi.noise_model(level)). A noise function is called as
noise(R, ids, phases, days, t_now) -> R, where days are the days each pair recorded together up to day t_now.
"""

import json
import time

import numpy as np
import pandas as pd

DATA = "../data"
DVP = 0.05                                      # fractional Vp change inside the BDT zone (+ = faster)
DOME = dict(top=(4.0, 6.5), base=(6.0, 9.0), width=2.5)   # km: depth at the crest and far away; Gaussian width
EDGE = 0.3                                      # km over which the change ramps in at the top and base


class MatrixImaging:
    """Inputs, settings and functions for P-wave matrix imaging of a dome-shaped BDT."""

    def __init__(self, data=DATA, scenarios=None, dvp=DVP, dome=None, edge=EDGE,
                 f_band=(1.0, 6.0), df=0.1,
                 image_dxy=0.25, image_dz=0.1, image_r=6.5, image_z=(0.5, 8.0),
                 scat_dxy=0.25, scat_dz=0.2, box=9.0, box_taper=2.0,
                 snr_min=3.0, top_tol=0.3, crest_frac=0.5, pick_z=(3.0, 8.0)):
        self.DATA = data
        self.SCENARIOS = {"dome": (0.0, 0.0)} if scenarios is None else dict(scenarios)   # name: BDT high (km)
        self.DVP, self.DOME, self.EDGE = dvp, dict(DOME, **(dome or {})), edge
        self.F = np.arange(f_band[0], f_band[1] + 1e-9, df)  # df small enough that two-way times do not wrap
        self.WF = np.hanning(len(self.F) + 2)[1:-1]          # source spectrum: Hann taper over the band
        self.SNR_MIN = snr_min                    # a column is imaged if its peak is SNR_MIN x the image median
        self.TOP_TOL = top_tol                    # km: a picked top within this of the truth counts as correct
        self.CREST_FRAC = crest_frac              # BDT high picked only from columns this bright vs the brightest
        self.PICK_Z = pick_z                      # km: the BDT top is picked only in this depth range (the image
                                                  # covers more; near-surface artefacts are not picked)
        self._load(image_dxy, image_dz, image_r, image_z, scat_dxy, scat_dz, box, box_taper)

    # ------------------------------------------------------------------------------------------------------------
    # Inputs and grids
    # ------------------------------------------------------------------------------------------------------------
    def _load(self, image_dxy, image_dz, image_r, image_z, scat_dxy, scat_dz, box, box_taper):
        d = self.DATA
        vm = np.load(f"{d}/velocity_model.npz")
        self.vx, self.vy, self.vz, self.v = vm["x"], vm["y"], vm["z"], vm["vp"]
        self.stations = pd.read_csv(f"{d}/stations.csv", index_col="station")
        self.D = json.load(open(f"{d}/designs.json"))
        tg = np.load(f"{d}/tt/grid.npz")
        self.tx, self.ty, self.tz = tx, ty, tz = tg["x"], tg["y"], tg["z"]
        # travel-time grid -> velocity-model grid indices (the tt grid is a crop of the model grid)
        self.mx, self.my, self.mz = (np.searchsorted(np.round(a, 6), np.round(b, 6))
                                     for a, b in ((self.vx, tx), (self.vy, ty), (self.vz, tz)))
        step = lambda km, g: max(1, int(round(km / (g[1] - g[0]))))

        # scatterer candidates: every node the BDT zone can reach, tapered towards the edge of the box
        zr = (self.DOME["top"][0] - self.EDGE, self.DOME["base"][1] + self.EDGE)
        si, sj, sk = self._nodes(box, zr, (step(scat_dxy, tx), step(scat_dxy, ty), step(scat_dz, tz)))
        self.SI, self.SJ, self.SK = (a.ravel() for a in np.meshgrid(si, sj, sk, indexing="ij"))
        edge = np.clip((box - np.maximum(np.abs(tx[self.SI]), np.abs(ty[self.SJ]))) / box_taper, 0, 1)
        self.SC_TAPER = 0.5 - 0.5 * np.cos(np.pi * edge)
        self.CELL = ((tx[1] - tx[0]) * step(scat_dxy, tx) * (ty[1] - ty[0]) * step(scat_dxy, ty)
                     * (tz[1] - tz[0]) * step(scat_dz, tz))

        # image grid
        self.ii, self.ij, self.ik = self._nodes(image_r, image_z, (step(image_dxy, tx), step(image_dxy, ty),
                                                                  step(image_dz, tz)))
        self.IX, self.IY, self.IZ = tx[self.ii], ty[self.ij], tz[self.ik]
        self.II, self.IJ, self.IK = (a.ravel() for a in np.meshgrid(self.ii, self.ij, self.ik, indexing="ij"))
        self.XX, self.YY = np.meshgrid(self.IX, self.IY, indexing="ij")

        # deployment
        self.PERMANENT = self.D["trillium"] + self.D.get("existing", [])
        self.DEPLOY_END = max(p["end"] for ph in self.D["designs"].values() for p in ph)
        self.P2 = self.D["p2"]
        self.AZR = np.radians(self.P2["az"])
        self.MT_R = self.D.get("mt", {}).get("radius", 4.0)
        self.BANK = {(i, j): s for s, i, j in zip(self.D["bank"]["ids"], self.D["bank"]["i"], self.D["bank"]["j"])}
        self.KEYS = np.array(list(self.BANK))
        self._tt, self._echo = {}, {}

    def __repr__(self):
        return (f"MatrixImaging: image {len(self.ii)} x {len(self.ij)} x {len(self.ik)} "
                f"({self.IZ[0]:.1f}-{self.IZ[-1]:.1f} km deep), {len(self.SI)} scatterer candidates, "
                f"{len(self.F)} frequencies, Vp {self.DVP:+.0%} in the BDT zone, scenarios {list(self.SCENARIOS)}")

    def _nodes(self, xr, zr, step):
        i = np.where(np.abs(self.tx) <= xr + 1e-9)[0][::step[0]]
        j = np.where(np.abs(self.ty) <= xr + 1e-9)[0][::step[1]]
        k = np.where((self.tz >= zr[0] - 1e-9) & (self.tz <= zr[1] + 1e-9))[0][::step[2]]
        return i, j, k

    def designs(self):
        """Design names in designs.json (phased ones are run adaptively)."""
        return list(self.D["designs"])

    def is_phased(self, design):
        return design in self.P2.get("kinds", {})

    # ------------------------------------------------------------------------------------------------------------
    # The Earth: dome-shaped BDT zone
    # ------------------------------------------------------------------------------------------------------------
    def centre(self, scenario):
        return np.asarray(self.SCENARIOS[scenario] if isinstance(scenario, str) else scenario, float)

    def bdt(self, scenario, x, y):
        """BDT zone top and base depth (km) at map points x, y (km)."""
        c = self.centre(scenario)
        w = 1 - np.exp(-((np.asarray(x) - c[0]) ** 2 + (np.asarray(y) - c[1]) ** 2) / self.DOME["width"] ** 2)
        (t0, t1), (b0, b1) = self.DOME["top"], self.DOME["base"]
        return t0 + (t1 - t0) * w, b0 + (b1 - b0) * w

    def true_top(self, scenario):
        """True BDT top on the image map grid."""
        return self.bdt(scenario, self.XX, self.YY)[0]

    def perturbation(self, scenario):
        """Scatterer nodes and strengths sigma = d(1/v^2) x cell volume x edge taper."""
        top, base = self.bdt(scenario, self.tx[self.SI], self.ty[self.SJ])
        zz = self.tz[self.SK]
        w = np.clip(np.minimum(zz - top, base - zz) / self.EDGE + 0.5, 0, 1)
        v0 = self.v[self.mx[self.SI], self.my[self.SJ], self.mz[self.SK]]
        sig = (1 / (v0 * (1 + self.DVP * w)) ** 2 - 1 / v0 ** 2) * self.CELL * self.SC_TAPER
        keep = np.abs(sig) > 1e-6 * np.abs(sig).max()
        return np.where(keep)[0], sig[keep].astype(np.float32)

    # ------------------------------------------------------------------------------------------------------------
    # Travel times, echoes and image
    # ------------------------------------------------------------------------------------------------------------
    def station_tt(self, s):
        """(travel times to all scatterer candidates, to all image nodes) for station s, from the 03 tables."""
        if s not in self._tt:
            t = np.load(f"{self.DATA}/tt/{s:03d}_p.npy", mmap_mode="r")
            self._tt[s] = (np.asarray(t[self.SI, self.SJ, self.SK], np.float32),
                           np.asarray(t[self.II, self.IJ, self.IK], np.float32))
        return self._tt[s]

    def spreading(self, s, x, y, z):
        r = np.sqrt((x - self.stations.x[s]) ** 2 + (y - self.stations.y[s]) ** 2 + z ** 2)
        return (1 / (4 * np.pi * np.maximum(r, 0.05))).astype(np.float32)

    def echoes(self, scenario, ids):
        """Noise-free reflection matrix R(f) between stations ids, (nF, N, N). Station pairs are computed once per
        scenario and kept, so designs that share stations (the Trilliums, phase 1) do not recompute them."""
        key = tuple(self.centre(scenario))
        if key not in self._echo:
            k, sig = self.perturbation(scenario)
            self._echo[key] = dict(k=k, sig=sig, pos={}, R=np.zeros((len(self.F), 64, 64), np.complex64),
                                   known=np.zeros((64, 64), bool))
        e = self._echo[key]
        for s in ids:
            e["pos"].setdefault(s, len(e["pos"]))
        if len(e["pos"]) > len(e["known"]):                    # grow the cache
            m = max(len(e["pos"]), 2 * len(e["known"]))
            R, kn = np.zeros((len(self.F), m, m), np.complex64), np.zeros((m, m), bool)
            c = len(e["known"])
            R[:, :c, :c], kn[:c, :c] = e["R"], e["known"]
            e["R"], e["known"] = R, kn
        sel = np.array([e["pos"][s] for s in ids])
        todo = ~e["known"][np.ix_(sel, sel)].all(axis=1)
        if todo.any():
            k = e["k"]
            T = np.stack([self.station_tt(s)[0][k] for s in ids]).astype(np.float64)
            A = np.stack([self.spreading(s, self.tx[self.SI[k]], self.ty[self.SJ[k]], self.tz[self.SK[k]])
                          for s in ids])
            G = A * np.exp(-2j * np.pi * self.F[0] * T)            # stepped frequency by frequency
            step = np.exp(-2j * np.pi * (self.F[1] - self.F[0]) * T)
            rows = sel[todo]
            for n in range(len(self.F)):
                blk = ((G[todo] * e["sig"]) @ G.T).astype(np.complex64)
                e["R"][n][np.ix_(rows, sel)] = blk
                e["R"][n][np.ix_(sel, rows)] = blk.T
                G *= step
            e["known"][np.ix_(rows, sel)] = e["known"][np.ix_(sel, rows)] = True
        return e["R"][:, sel][:, :, sel]

    def image(self, ids, R, pair_mask):
        """Confocal image I(x) = |sum_f W(f) sum_ij G*(x, s_i) R_ij G*(s_j, x)| on the image grid."""
        T = np.stack([self.station_tt(s)[1] for s in ids]).T
        A = np.stack([self.spreading(s, self.tx[self.II], self.ty[self.IJ], self.tz[self.IK]) for s in ids]).T
        I = np.zeros(len(self.II), np.complex128)
        for n, f in enumerate(self.F):
            Gc = A * np.exp(2j * np.pi * f * T)
            I += self.WF[n] * np.einsum("pi,ij,pj->p", Gc, R[n] * pair_mask, Gc, optimize=True)
        return np.abs(I).reshape(len(self.ii), len(self.ij), len(self.ik)).astype(np.float32)

    # ------------------------------------------------------------------------------------------------------------
    # Who records with whom
    # ------------------------------------------------------------------------------------------------------------
    def windows(self, phases, t_now=None):
        """Recording windows (days) up to t_now (default: the end of the deployment): permanent stations
        throughout, nodes during their phase."""
        t_now = self.DEPLOY_END if t_now is None else t_now
        w = {s: [(0, t_now)] for s in self.PERMANENT}
        for p in phases:
            for s in p["stations"]:
                w.setdefault(s, [])
                if min(p["end"], t_now) > p["start"]:
                    w[s].append((p["start"], min(p["end"], t_now)))
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

    def image_phases(self, scenario, phases, noise=None, t_now=None):
        """Image from the data recorded up to day t_now (default: the end of the deployment) by every station in
        the phases plus the permanent ones. Returns (image, station ids, days)."""
        t_now = self.DEPLOY_END if t_now is None else t_now
        ids = sorted(set(self.PERMANENT) | {s for p in phases for s in p["stations"]})
        days = self.overlap_days(ids, self.windows(phases, t_now))
        R = self.echoes(scenario, ids)
        if noise is not None:
            R = noise(R, ids, phases, days, t_now)
        return self.image(ids, R, days > 0), ids, days

    def signal(self, scenario, design):
        """rms echo (noise-free |R|) over the station pairs of a static design: the reference for noise levels."""
        ph = self.D["designs"][design]
        ids = sorted(set(self.PERMANENT) | {s for p in ph for s in p["stations"]})
        mask = self.overlap_days(ids, self.windows(ph)) > 0
        return float(np.sqrt(np.mean(np.abs(self.echoes(scenario, ids)[:, mask]) ** 2)))

    def noise_model(self, level, days_ref=30, real=0, block_days=7):
        """Correlation noise with std level x sqrt(days_ref / days) for a pair that recorded together for `days`.

        Each pair's noise is the stack of its weekly records, drawn from (real, pair) alone: a phase 1 pick and the
        final image see the same noise on the data they share, and so do designs that share stations.
        """
        edges = np.arange(int(np.ceil(self.DEPLOY_END / block_days)) + 1) * block_days
        lo, hi = edges[:-1], edges[1:]
        nf = len(self.F)

        def shared(w1, w2):
            return [(max(a1, a2), min(b1, b2)) for a1, b1 in w1 for a2, b2 in w2 if min(b1, b2) > max(a1, a2)]

        def noise(R, ids, phases, days, t_now):
            win = self.windows(phases, t_now)
            acc = np.zeros(R.shape, np.complex64)
            for p, q in zip(*np.nonzero(np.triu(days > 0, 1))):
                d = sum(np.clip(np.minimum(hi, b) - np.maximum(lo, a), 0, None)
                        for a, b in shared(win[ids[p]], win[ids[q]]))
                nz = np.nonzero(d)[0]
                Z = np.random.default_rng([real, *sorted((int(ids[p]), int(ids[q])))]).standard_normal(
                    (nz[-1] + 1, 2, nf), dtype=np.float32)            # the pair's weekly record, first weeks only
                acc[:, p, q] = ((Z[nz, 0] + 1j * Z[nz, 1]) * np.sqrt(d[nz] / 2)[:, None]).sum(axis=0)
            acc += acc.transpose(0, 2, 1)
            std = level * np.sqrt(days_ref) / np.maximum(days, 1e-9) * (days > 0)
            return (R + acc * std).astype(np.complex64)

        return noise

    # ------------------------------------------------------------------------------------------------------------
    # Picks and scores
    # ------------------------------------------------------------------------------------------------------------
    def pick_top(self, img, min_peak=None):
        """Picked BDT top per image column: shallowest local maximum above half the column peak, within PICK_Z.
        NaN if the column peak is below the detection floor (SNR_MIN x the image median, or min_peak)."""
        kz = (self.IZ >= self.PICK_Z[0] - 1e-9) & (self.IZ <= self.PICK_Z[1] + 1e-9)
        sub, z = img[..., kz], self.IZ[kz]
        peak = sub.max(axis=2)
        imaged = peak >= (self.SNR_MIN * np.median(img) if min_peak is None else min_peak)
        a = sub / np.maximum(peak[..., None], 1e-30)
        local = np.zeros_like(a, bool)
        local[..., 1:-1] = (a[..., 1:-1] >= a[..., :-2]) & (a[..., 1:-1] >= a[..., 2:]) & (a[..., 1:-1] >= 0.5)
        first = np.argmax(local, axis=2)
        return np.where(imaged & local.any(axis=2), z[first], np.nan)

    def pick_target(self, img, tol=0.1):
        """BDT high from an image: centroid of the bright columns (peak >= CREST_FRAC x the image max) whose top
        is within tol km of the shallowest one. None if nothing is picked."""
        top = self.pick_top(img, self.CREST_FRAC * img.max())
        if not np.isfinite(top).any():
            return None
        hi = np.isfinite(top) & (top <= np.nanmin(top) + tol)
        return np.array([self.XX[hi].mean(), self.YY[hi].mean()])

    def score(self, img, scenario):
        """Fraction of the MT target area whose BDT top is imaged within TOP_TOL, the imaged depth at the true
        crest, and how far the picked BDT high is from the true one."""
        c = self.centre(scenario)
        true, est = self.true_top(scenario), self.pick_top(img)
        area = np.hypot(self.XX, self.YY) <= self.MT_R
        good = area & (np.abs(est - true) <= self.TOP_TOL)
        i, j = np.argmin(np.abs(self.IX - c[0])), np.argmin(np.abs(self.IY - c[1]))
        kz = (self.IZ >= self.PICK_Z[0] - 1e-9) & (self.IZ <= self.PICK_Z[1] + 1e-9)
        col = img[i, j, kz]
        target = self.pick_target(img)
        return dict(area_coverage=float(good.sum() / area.sum()),
                    area_imaged=float((area & np.isfinite(est)).sum() / area.sum()),
                    crest_true_km=float(true[i, j]),
                    crest_imaged_km=float(self.IZ[kz][np.argmax(col)]),
                    target_error_km=float(np.hypot(*(target - c))) if target is not None else np.nan)

    # ------------------------------------------------------------------------------------------------------------
    # Phase 2 around a target, with 02's rules
    # ------------------------------------------------------------------------------------------------------------
    def _block_index(self, centre, size):
        P2 = self.P2
        d = np.asarray(centre, float)
        if np.hypot(*d) > P2["search_r"]:
            d *= P2["search_r"] / np.hypot(*d)
        a = self.AZR
        along, across = d[0] * np.sin(a) + d[1] * np.cos(a), -d[0] * np.cos(a) + d[1] * np.sin(a)
        out = []
        for c, n in zip((along, across), size):
            frac = ((n - 1) / 2) % 1
            base = np.round(c / P2["spacing"] - frac)
            out.append(np.rint(base + frac + np.arange(n) - (n - 1) / 2).astype(int))
        return out

    def phase2_sites(self, centre, kind):
        """Bank (i, j) sites of each phase 2 patch around centre, in deployment order (as 02's phase2)."""
        (pa, pc), (na, nc) = self.P2["patch"], self.P2["patches"]
        I, J = self._block_index(centre, (pa * na, pc * nc))
        if kind == "rolling":
            order = [(a, b) for a in range(na) for b in (range(nc) if a % 2 == 0 else reversed(range(nc)))]
            return [[(i, j) for i in I[pa * a:pa * (a + 1)] for j in J[pc * b:pc * (b + 1)]] for a, b in order]
        if kind == "interleaved":
            K = na * nc
            cls = (lambda a, b: (a + b) % K) if nc == 1 else (lambda a, b: (a % na) * nc + b % nc)
            return [[(I[a], J[b]) for a in range(len(I)) for b in range(len(J)) if cls(a, b) == q]
                    for q in range(K)]
        raise ValueError(kind)

    def bank_ids(self, patches):
        """Station ids of bank sites; a site in water is replaced by the nearest free dry bank site (as 02's ids)."""
        used = {s for p in patches for s in p}
        out = []
        for p in patches:
            out.append([])
            for site in p:
                if site not in self.BANK:
                    order = np.argsort(np.hypot(self.KEYS[:, 0] - site[0], self.KEYS[:, 1] - site[1]), kind="stable")
                    site = next(k for k in map(lambda r: (int(r[0]), int(r[1])), self.KEYS[order]) if k not in used)
                    used.add(site)
                out[-1].append(self.BANK[site])
        return out

    # ------------------------------------------------------------------------------------------------------------
    # One design, one scenario
    # ------------------------------------------------------------------------------------------------------------
    def run_design(self, design, scenario, noise=None):
        """Image one design. Phased designs: image the data recorded by the end of phase 1, pick the BDT high, lay
        phase 2 out around it.

        Returns dict(image, phases, ids, score, seconds) and, for phased designs, phase1 (its image) and target.
        """
        t0 = time.time()
        ph = self.D["designs"][design]
        out = {}
        if self.is_phased(design):
            p1 = ph[0]
            img1, _, _ = self.image_phases(scenario, [p1], noise, t_now=p1["end"])   # data recorded by then
            target = self.pick_target(img1)
            out.update(phase1=img1, target=target if target is not None else np.array(self.P2["centre0"], float))
            L = self.P2["patch_days"]
            ph = [p1] + [dict(start=p1["end"] + q * L, end=p1["end"] + (q + 1) * L, stations=st)
                         for q, st in enumerate(self.bank_ids(self.phase2_sites(out["target"],
                                                                                self.P2["kinds"][design])))]
        img, ids, days = self.image_phases(scenario, ph, noise)
        out.update(image=img, phases=ph, ids=ids, n_pairs=int((days > 0).sum() // 2),
                   score=self.score(img, scenario), seconds=time.time() - t0)
        return out

    # ------------------------------------------------------------------------------------------------------------
    # Plot helpers
    # ------------------------------------------------------------------------------------------------------------
    def line(self, scenario, half=6.0, n=241):
        """Points along the node axis (AZ) through the scenario's BDT high: (s, x, y)."""
        c = self.centre(scenario)
        s = np.linspace(-half, half, n)
        return s, c[0] + s * np.sin(self.AZR), c[1] + s * np.cos(self.AZR)

    def section(self, img, scenario, half=6.0):
        """Image along the line, normalised to the image max (NaN outside the image). Returns (s, z, section)."""
        s, x, y = self.line(scenario, half)
        inside = (x >= self.IX[0]) & (x <= self.IX[-1]) & (y >= self.IY[0]) & (y <= self.IY[-1])
        jx = np.clip(np.searchsorted(self.IX, x), 0, len(self.IX) - 1)
        jy = np.clip(np.searchsorted(self.IY, y), 0, len(self.IY) - 1)
        return s, self.IZ, np.where(inside[:, None], img[jx, jy, :] / img.max(), np.nan)

    def velocity_section(self, scenario, half=6.0):
        """Vp along the line with the BDT zone applied. Returns (s, z, vp, top, base)."""
        s, x, y = self.line(scenario, half)
        kx = np.clip(np.round((x - self.vx[0]) / (self.vx[1] - self.vx[0])).astype(int), 0, len(self.vx) - 1)
        ky = np.clip(np.round((y - self.vy[0]) / (self.vy[1] - self.vy[0])).astype(int), 0, len(self.vy) - 1)
        top, base = self.bdt(scenario, x, y)
        w = np.clip(np.minimum(self.vz[None, :] - top[:, None], base[:, None] - self.vz[None, :]) / self.EDGE + 0.5,
                    0, 1)
        return s, self.vz, self.v[kx, ky, :] * (1 + self.DVP * w), top, base
