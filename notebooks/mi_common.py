"""Shared setup for the matrix imaging notebooks (05, 06, 07): inputs, settings and functions."""

# --- Imports ---
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import map_coordinates

# --- Stations and designs ---
stations = pd.read_csv("../data/stations.csv", index_col="station")
D = json.load(open("../data/designs.json"))
designs, bank, p2 = D["designs"], D["bank"], D["p2"]

# --- Travel-time tables ---
g = np.load("../data/tt/grid.npz")
dg = np.array([g["x"][1] - g["x"][0], g["y"][1] - g["y"][0], g["z"][1] - g["z"][0]])   # table cell size
g0 = np.array([g["x"][0], g["y"][0], g["z"][0]])                                        # table corner
t1 = np.load("../data/tt/tt1d.npz")
tt1 = RectBivariateSpline(t1["r"], t1["z"], t1["t_p"])

# --- Scenarios ---
S = json.load(open("../data/scenarios.json"))
sc = np.load("../data/scenarios.npz")
mx, my = sc["x"], sc["y"]
MX, MY = np.meshgrid(mx, my, indexing="ij")

# --- Settings ---
FREQS = np.linspace(2, 6, 30)          # frequencies used (Hz)
REFS = ["constant", "1d", "3d"]        # focusing models
REF_DECIDE = "1d"                      # model used to make the phase 1 decision (the realistic one)
C0 = 4.5                               # constant focusing velocity (km/s)
DROP, ROUGH = 0.05, 0.3                # Vp drop into the BDT (sets the reflection) and patchiness
SCAT_R = 7.0                           # scatterers on the true BDT top within this distance of the conductor (km)
STATIC = ["grid", "wide grid", "spiral"]

P1_H = np.arange(-7, 7.01, 0.5)        # phase 1 image: map positions (km)
P1_Z = np.arange(2, 9.01, 0.25)        # phase 1 image: depths (km)

R_EVAL = 3.0                           # score the BDT top within this distance of the conductor (km)
BRIGHT, FOUND_DZ = 0.3, 0.5            # a column counts if bright enough and within 0.5 km of the true top

rng = np.random.default_rng(4)
r0 = -DROP / (2 - DROP)                # reflection off fast -> slow (same density)


# --- Travel times and focusing ---
def tt3d(st, pts):
    """3D travel times (s) from stations st to points, from the tables."""
    idx = ((pts - g0) / dg).T
    return np.array([map_coordinates(np.load(f"../data/tt/{k:03d}_p.npy").astype(np.float32), idx, order=1)
                     for k in st])


def dist(st, pts):
    """Straight-line distance (km) from stations st to points (at least 0.1 km)."""
    return np.maximum(0.1, np.sqrt((stations.x.values[st, None] - pts[:, 0]) ** 2 +
                                   (stations.y.values[st, None] - pts[:, 1]) ** 2 + pts[:, 2] ** 2))


def tt_ref(st, pts, ref):
    """Travel times from a focusing model: '3d', '1d' or 'constant'."""
    if ref == "3d":
        return tt3d(st, pts)
    if ref == "constant":
        return dist(st, pts) / C0
    rr = np.hypot(stations.x.values[st, None] - pts[:, 0], stations.y.values[st, None] - pts[:, 1])
    return tt1.ev(rr, np.broadcast_to(pts[:, 2], rr.shape))


def reflection_matrix(st, scat, sigma):
    """R_ij(f): echoes from every scatterer, through the true 3D model."""
    T, A = tt3d(st, scat), 1 / (4 * np.pi * dist(st, scat))
    return np.array([(A * np.exp(-2j * np.pi * f * T) * sigma) @ (A * np.exp(-2j * np.pi * f * T)).T
                     for f in FREQS])


def image(R, st, pts, ref):
    """Confocal image: undo the focusing model's delay for every station pair and point, then add up."""
    T, A = tt_ref(st, pts, ref), 1 / (4 * np.pi * dist(st, pts))
    I = np.zeros(len(pts), complex)
    for Rf, f in zip(R, FREQS):
        G = A * np.exp(-2j * np.pi * f * T)
        I += np.einsum("ip,ij,jp->p", G.conj(), Rf, G.conj())
    return np.abs(I)


# --- Phase 2 from the station bank ---
ua = np.array([np.sin(np.radians(p2["az"])), np.cos(np.radians(p2["az"]))])   # along AZ
un = np.array([-np.cos(np.radians(p2["az"])), np.sin(np.radians(p2["az"]))])  # across AZ
bid, bi, bj = np.array(bank["ids"]), np.array(bank["i"]), np.array(bank["j"])
ALONG, ACROSS = np.arange(15) - 7, np.arange(8) - 4          # rolling grid, in bank steps


def phase2_patches(centre):
    """Station ids of each phase 2 patch, centred as close as the bank allows to centre."""
    d = np.asarray(centre) - np.asarray(p2["centre0"])
    ci, cj = round(d @ ua / p2["spacing"]), round(d @ un / p2["spacing"])
    return [bid[np.isin(bi, ci + ALONG[ia * 5:(ia + 1) * 5]) & np.isin(bj, cj + ACROSS[ja * 4:(ja + 1) * 4])].tolist()
            for ia, ja in p2["order"]]


# --- Image points and scoring ---
V1x, V1y, V1z = np.meshgrid(P1_H, P1_H, P1_Z, indexing="ij")
p1_pts = np.column_stack([V1x.ravel(), V1y.ravel(), V1z.ravel()])


def columns(centre):
    """Vertical columns on a 0.5 km map grid around a centre, for scoring."""
    h = np.arange(-R_EVAL, R_EVAL + 0.01, 0.5)
    cx, cy = [a.ravel() for a in np.meshgrid(centre[0] + h, centre[1] + h, indexing="ij")]
    zc = np.arange(2, 9.01, 0.2)
    pts = np.column_stack([np.repeat(cx, len(zc)), np.repeat(cy, len(zc)), np.tile(zc, len(cx))])
    return cx, cy, zc, pts


def top_at(top, x, y):
    """True BDT top depth at map points (nearest model cell)."""
    i = np.clip(np.round((x - mx[0]) / (mx[1] - mx[0])).astype(int), 0, len(mx) - 1)
    j = np.clip(np.round((y - my[0]) / (my[1] - my[0])).astype(int), 0, len(my) - 1)
    return top[i, j]


def coverage(I, cx, cy, zc, true_top, centre):
    """Fraction of the true top (within R_EVAL of the conductor) found: bright column, within FOUND_DZ."""
    I = I.reshape(len(cx), len(zc))
    I = I / I.max()
    use = ~np.isnan(true_top) & (np.hypot(cx - centre[0], cy - centre[1]) <= R_EVAL)
    found = (I.max(axis=1) >= BRIGHT) & (np.abs(zc[I.argmax(axis=1)] - true_top) <= FOUND_DZ)
    return np.mean(found[use])