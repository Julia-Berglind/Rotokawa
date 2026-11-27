"""
geonet_matrix_imaging.py
Download one hour of vertical-component GeoNet data for the five stations around Rotokawa, save it, and run the
matrix imaging pipeline with the Rotokawa velocity model.

Run from notebooks/ (or anywhere: paths are relative to the repo's data/ folder):
    python geonet_matrix_imaging.py

Inputs:  data/site.json, data/velocity_model.npz (notebook 01)
Outputs: data/geonet_matrix_imaging/
           waveforms/NZ.<STA>..<CHA>.sac    response-removed velocity, 20 Hz, with stla/stlo/stel set
           inventory.xml                    station metadata with responses
           mi_results.npz, mi_correlations.png, mi_images.png, mi_phase_velocity.png
"""
from obspy import UTCDateTime, read, read_inventory
from obspy.clients.fdsn import Client

import matrix_imaging_pipeline as mip

STATIONS = "WPRZ,ARAZ,THQ2,ALRZ,THSZ"
CHANNELS = "HHZ,EHZ"
T0 = UTCDateTime("2026-11-20T14:00:00")
T1 = T0 + 3600
FS = 20.0                                    # resample to this (Hz)
PRE_FILT = (0.05, 0.1, 8, 9.5)               # Hz, for response removal
REDOWNLOAD = False                           # True: fetch again even if waveforms/ already has files

OUT = mip.OUTDIR
WAV = OUT / "waveforms"
WAV.mkdir(parents=True, exist_ok=True)
INV = OUT / "inventory.xml"

if REDOWNLOAD or not list(WAV.glob("NZ.*.sac")):
    c = Client("GEONET")
    inv = c.get_stations(network="NZ", station=STATIONS, level="channel", starttime=T0, endtime=T1)
    for net in inv:
        for sta in net:
            chans = sorted({f"{ch.location_code}.{ch.code} ({ch.sample_rate:g} Hz)" for ch in sta})
            print(sta.code, f"{sta.latitude:.4f} {sta.longitude:.4f}", chans)

    inv = c.get_stations(network="NZ", station=STATIONS, channel=CHANNELS, level="response",
                         starttime=T0, endtime=T1)
    inv.write(str(INV), format="STATIONXML")
    st = c.get_waveforms("NZ", STATIONS, "*", CHANNELS, T0 - 60, T1 + 60)
    print(st)                                # check which stations came back
    st.merge(fill_value=0)
    st.remove_response(inventory=inv, output="VEL", pre_filt=PRE_FILT)
    st.resample(FS).trim(T0, T1)
    for tr in st:
        co = inv.get_coordinates(tr.id, T0)
        tr.stats.sac = {"stla": co["latitude"], "stlo": co["longitude"], "stel": co["elevation"]}
        tr.write(str(WAV / f"{tr.id}.sac"), format="SAC")

# read every saved station, HHZ and EHZ alike (the pipeline keeps one channel per station, HHZ first)
st = read(str(WAV / "NZ.*.sac"))
print(st)
out = mip.run(st, vmodel=mip.VMODEL, outdir=OUT)
