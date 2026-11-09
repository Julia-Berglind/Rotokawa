# Matrix Imaging Forward Modelling Notes
My notes on forward modelling with the matrix imaging approach.

A seismometer records ground motion against time. If a pulse leave a source and reaches a station along several paths, the record is the pulse repeated. Each copy is deplayed by that path's travel time and scaled by how much energy survives the trip. We can express this as:

$$
u(t) = \sum_{\text{paths}} a_p \, w(t - \tau_p)
$$
 
- $w(t)$: the pulse (wavelet)
- $\tau_p$: travel time of path $p$
- $a_p$: amplitude of path $p$

Everything that follows is about working out the $\tau_p$ and $a_p$ for echoes off of the BDT, then using them in reverse to find out where the echoes came from.

First, travel time. In a rock of constant velocity, $v$, a wave travels in a straight line:

$$
\tau = \frac{r}{v}
$$

where $r$ is the distance. In real rock, the velocity changes from place to place, and waves bend. The travel time, $T(X)$, from a source then obeys the **eikonal equation**:

$$
|\Delta T| = \frac{1}{v(x)}
$$

This tells us that the travel time grows fastest in the direction the wave travels, at a rate of one over the local speed. '03_travel_times' gives the first-arrival from a station to every point in the model. For now, we will keep $v$ contsant so $\tau = r/v$

In terms of the amplitude, we can apply the concept of geometric spreading. A wave that leaves a point spreads over a sphere of area $4\pi r^2$. Energy is conserved, so unit area falls as $1/r^2$ and amplitude (square root of energy) falls as $1/r$. With the standard normalisation:

$$
a = \frac{1}{4 \pi r}
$$

Putting these two concepts together, the response at $\mathbf{x}$ to a perfectly sharp pulse $\delta(t)$ fired at $s$ is **Green's function**:

$$
G(\mathbf{x}, \mathbf{s}, t) = \frac{\delta(t - r/v)}{4\pi r}, \qquad r = |\mathbf{x} - \mathbf{s}|
$$

We then shift this into the frequency domain using a fourier transform to see how much of each frequency it contains. Due to the shift property, delaying a signal by $\tau$ mutliples its transform by a phase factor:

$$w(t - \tau) \;\longrightarrow\; W(f)\, e^{-i 2\pi f \tau}$$

This makes a delay in the time domain become multipliation in the frequency domain. Thus, Green's function in the frequency domain becomes:


$$
G(\mathbf{x}, \mathbf{s}, f) = \frac{e^{-i 2\pi f r / v}}{4\pi r}
$$

This is important: **in the frequency domain, delays are rotations!**

If we consider a normal scatterer, a small spot where the rock's properties change suddenly, when a wave hits it, the spot acts as a tiny new source and then re-radiates some of the energy in all direction. Its strength, $\sigma$ says how much. 

For a wave sent from station $i$, scattered at point $\mathbf{k}$, and recorded at station $j$:
 
1. Station $i$ to point $\mathbf{k}$: Green's function $G(\mathbf{k}, \mathbf{s}_i)$
2. Re-radiated with strength $\sigma_k$
3. Point $\mathbf{k}$ to station $j$: Green's function $G(\mathbf{s}_j, \mathbf{k})$

In the ferquency domain, these steps multiply in sequence such that:

$$
u_{ij}(f) = G(\mathbf{s}_j, \mathbf{k}, f)\; \sigma_k \; G(\mathbf{k}, \mathbf{s}_i, f)
$$
 
Writing the Green's functions out:
 
$$
u_{ij}(f) = \frac{\sigma_k}{(4\pi r_{ik})(4\pi r_{kj})}\; e^{-i 2\pi f (r_{ik} + r_{kj})/v}
$$

The phases add, so the echo arrives at a total travel time $r_{ik} + r_{kj})/v$ and the amplitude is weakened by spreading on both legs.

Now we apply the **born approximation** which assumes the wave scatteres once and that the scatterers are weak such that they send a little energy back but do not change the travel times of waves passing through. This approximation is valid for small contrasts, such as a few percent change in velocity, which is what a BDT low-velocity zone is likely to be.

With one scatterer, the math is fairly clean, but in the ground, we have many scatterers, and waves add. This means each scatterer sends it's own echo and the record we get is just the sum:

$$
R_{ij}(f) = \sum_k G(\mathbf{s}_j, \mathbf{k}, f)\; \sigma_k \; G(\mathbf{k}, \mathbf{s}_i, f)
$$

$R_{ij}(f)$ is the complete response at station $j$ to a source at station $i$, at frequency $f$.

At a reflecting surface, like the top of the BDT, we can model it as points on that surface each with a strength reflection coefficient. If the points are much closer together than a wavelegnth, then waves can't tell them apart and 'see' a continuous surface. Their echoes interfere to make a proper reflection.

So, we come to this idea of a relfection matrix. We collect $R_{ij}$ for every station pair and put it into a matrix, $R(f)$ with $N$ stations, it is $N \times N$, one per frequency.

At this point, we could write $G$ for the matrix of Green's function for every station (rows) to every scatterer (columns) and $\Sigma$ for the diagonal matrix of the $\sigma_k$, then the sum over $k$ is a matrix product:

$$
\mathbf{R}(f) = \mathbf{G}(f)\ \boldsymbol{\Sigma}\ \mathbf{G}(f)^{\mathsf{T}}
$$

$R$ is symmetric because swapping source and reiever reverses the path without changing its length.

In real life, and at Rotokawa, there are no controlled sources. Instead, ambient noise from ocean, wind, traffic, etc arrives from all different directions. Cross-correlating the noise recorded at station $i$ and $j$, over a long time, converges to the Green's function between them, as if $i$ were a source. This is the concept of seismic noise interferometry.

In my code, I skip this step and calculate R directly, which is the result of perfectly converged correlctions. '06_noise' adds leftover noise back.

Now the idea is that we have $R$ and we want to know where the scatterers are. So, we pick a point $x$ underground and ask, "If a scatterer were at $x$, it's echo pair $(i, j)$ would carry the phase $G(x, s_i) G(s_j, x)$. If I undo that phase on $R_{ij}$, do the echoes for all pairs line up?"

Undo the predicted phase for every pair, add over all pairs and frequencies, and take the length of the result:
 
$$
I(\mathbf{x}) = \left| \sum_f \sum_{i,j} G^* (\mathbf{x}, \mathbf{s}_i, f)\; R_{ij}(f)\; G^* (\mathbf{s}_j, \mathbf{x}, f) \right|
$$

This is called the confocal image.

## Controls on image quality
The first control is frequency. At one frequency, a wrong point can line up by accident when its leftover delay is a whole number of periods. This gives false bright spots. However, a leftover delay $\delta$ rotates by $2 \pi f \delta$ which is a different angle at every frequency. Only the true point will line up at all of them. So, summing many frequencies keeps the true point and washes out the stripes.

The image also uses predicted phases. If the velocity used to predict is wrong, the leftover delays are not zero at the true point, so the image shifts and blurs. A fast focusing velocity puts the reflector too deep.
