# Klimontovich Equation — 1D–1V Numerical Solver

## What is the Klimontovich Equation?

The Klimontovich equation is the **exact kinetic equation** for a classical plasma. It tells us about the time evolution of the phase-space distribution function $f(\mathbf{x}, \mathbf{v}, t)$, defined as a sum of delta functions over all $N$ particles:

$$f(\mathbf{x}, \mathbf{v}, t) = \sum_{a=1}^{N} \delta\!\bigl(\mathbf{x} - \mathbf{x}_a(t)\bigr)\,\delta\!\bigl(\mathbf{v} - \mathbf{v}_a(t)\bigr)$$

It was formulated by **A. A. Klimontovich** (1961) as the well considered as starting point for plasma kinetic theory — exact for any finite $N$, from which the Vlasov equation, Vlasov–Poisson, Vlasov–Maxwell, and Landau/Balescu–Lenard collision operators are recovered in appropriate limits.

## The Equation (1D–1V Form)

$$\frac{\partial f}{\partial t} + v\,\frac{\partial f}{\partial x} + \frac{qE}{m}\,\frac{\partial f}{\partial v} = 0$$

| Symbol | Meaning | Value (this code) |
|--------|---------|-------------------|
| $q$ | Particle charge | $+1.602176634 \times 10^{-19}$ C (proton) |
| $m$ | Particle mass | $1.676221923 \times 10^{-27}$ kg (proton) |
| $E$ | Electric field | User-specified, constant (V/m) |
| $x$ | Position coordinate | $[0, L]$, periodic |
| $v$ | Velocity coordinate | $[-v_{\max},\, +v_{\max}]$ |

The 3D form includes $\mathbf{v}\cdot\nabla_{\mathbf{x}} f$ and $\frac{q}{m}(\mathbf{E}+\mathbf{v}\times\mathbf{B})\cdot\nabla_{\mathbf{v}} f$; this code solves the 1D–1V reduction with a prescribed (non-self-consistent) field.

## Numerical Scheme

| Term | Method | Order |
|------|--------|-------|
| $v\,\partial f/\partial x$ | Upwind (direction depends on sign of $v$) | 1st |
| $\frac{qE}{m}\,\partial f/\partial v$ | Central difference | 2nd |
| $\partial f/\partial t$ | Forward Euler | 1st |

Time step: $\Delta t = 0.5 \times \dfrac{\Delta x}{2\,v_{\max}}$ (CFL-limited with safety factor).

Boundary conditions: **periodic** in $x$, **zero-gradient** in $v$.

Initial condition: Gaussian bump $f \propto \exp\!\bigl[-(x-L/2)^2/\sigma_x^2\bigr]\,\exp\!\bigl[-v^2/v_{th}^2\bigr]$.

## Compilation and Execution

**Prerequisites:** GNU Fortran compiler (`gfortran`) version 4.0 or later.

**suggested input**
Nx:        256
Nv:        256
L:         1.0
v_max:     1.0e5
E_ext:     1.0e3
nsteps:    1000
sigma_x:   0.05
v_th:      1.0e4   

**OUTPUT(varies depends on input values)**
 ============================================================
   KLIMONTOVICH EQUATION  (1D-1V, prescribed E-field)
   df/dt + v*df/dx + (qE/m)*df/dv = 0
 ============================================================
 
 Enter grid size in x (Nx):
256
 Enter grid size in v (Nv):
1024
 Enter domain length L (m):
1.0
 Enter max velocity v_max (m/s):
1.0e6
 Enter external E-field (V/m):
1.0e3
 Enter number of time steps:
1000
 Enter bump width sigma_x (m):
0.05
 Enter thermal velocity v_th (m/s):
1.0e4
 
 dx =    3.9062500000000000E-003  m
 dv =    1953.1250000000000       m/s
 dt =    9.7656249999999996E-010  s  (CFL-limited)
 
 Step            1   N =    1570.7963267948871        <v> =    93.342390865549277        U =    39413093604.155350     
 Step          100   N =    1570.7963267948860        <v> =    9334.2390865546022        U =    121334336092.87811     
 Step          200   N =    1570.7963267948801        <v> =    18668.478173109201        U =    340259136732.44696     
 Step          300   N =    1570.7963267948869        <v> =    28002.717259663830        U =    696044310088.58191     
 Step          400   N =    1570.7963267948805        <v> =    37336.956346218292        U =    1188689856161.2883     
 Step          500   N =    1570.7963267948785        <v> =    46671.195432772962        U =    1818195774950.5522     
 Step          600   N =    1570.7963267948726        <v> =    56005.434519327508        U =    2584562066456.3916     
 Step          700   N =    1570.7963267948862        <v> =    65339.673605882119        U =    3487788730678.7842     
 Step          800   N =    1570.7963267948812        <v> =    74673.912692436745        U =    4527875767617.7920     
 Step          900   N =    1570.7963267948744        <v> =    84008.151778991407        U =    5704823177273.2754     
 Step         1000   N =    1570.7963267948892        <v> =    93342.390865544541        U =    7018630959645.3721     
 
 ============================================================
   FINAL DIAGNOSTICS
 ============================================================
 Total particles (integral of f):    1570.7963267948892     
 Mean velocity <v> (m/s):             93342.390865544541     
 Kinetic energy U (arb. units):       7018630959645.3721     
 
 Note: <v> should shift by (qE/m)*t over the simulation.
 Expected shift = (q*E/m)*t_total =    93342.390865545603       m/s


```bash
gfortran klimontovich.f90 -o klimontovich
./klimontovich   


