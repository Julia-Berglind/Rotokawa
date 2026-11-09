# Matrix Imaging Forward Modelling Notes

A seismic record is a sum of delayed, scaled copies, which can be mathematically represented as:

$$
u(t) = \sum_{\text{paths}} a_p \, w(t - \tau_p)
$$
 
- $w(t)$: the pulse (wavelet)
- $\tau_p$: travel time of path $p$
- $a_p$: amplitude of path $p$
