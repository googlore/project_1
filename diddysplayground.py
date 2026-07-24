import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.fft import fft, fftfreq

# ====================================================
# Advanced Computational Physics Simulation
# Nonlinear Driven Duffing Oscillator
# ====================================================

# Physical parameters
m = 1.0
k = 1.0
beta = 0.25      # Nonlinear spring coefficient
gamma = 0.15     # Damping
F = 0.35         # Driving amplitude
omega = 1.2      # Driving frequency

dt = 0.01
T = 200
N = int(T / dt)

# Initial conditions
x = np.zeros(N)
v = np.zeros(N)
t = np.linspace(0, T, N)

x[0] = 1.0
v[0] = 0.0


# ====================================================
# Differential equations
# ====================================================

def acceleration(x, v, time):
    return (
        -k*x
        -beta*x**3
        -gamma*v
        +F*np.cos(omega*time)
    ) / m


# ====================================================
# RK4 Integrator
# ====================================================

for i in range(N-1):

    k1x = v[i]
    k1v = acceleration(x[i], v[i], t[i])

    k2x = v[i] + 0.5*dt*k1v
    k2v = acceleration(
        x[i]+0.5*dt*k1x,
        v[i]+0.5*dt*k1v,
        t[i]+0.5*dt
    )

    k3x = v[i] + 0.5*dt*k2v
    k3v = acceleration(
        x[i]+0.5*dt*k2x,
        v[i]+0.5*dt*k2v,
        t[i]+0.5*dt
    )

    k4x = v[i] + dt*k3v
    k4v = acceleration(
        x[i]+dt*k3x,
        v[i]+dt*k3v,
        t[i]+dt
    )

    x[i+1] = x[i] + dt*(k1x+2*k2x+2*k3x+k4x)/6
    v[i+1] = v[i] + dt*(k1v+2*k2v+2*k3v+k4v)/6


# ====================================================
# Energy
# ====================================================

KE = 0.5*m*v**2
PE = 0.5*k*x**2 + 0.25*beta*x**4
Energy = KE + PE


# ====================================================
# FFT Analysis
# ====================================================

yf = np.abs(fft(x))
xf = fftfreq(N, dt)

positive = xf > 0


# ====================================================
# Poincare Section
# ====================================================

period = 2*np.pi/omega

indices = []

for i in range(N):
    if abs((t[i] % period)) < dt:
        indices.append(i)

px = x[indices]
pv = v[indices]


# ====================================================
# Plot Layout
# ====================================================

fig = plt.figure(figsize=(16,10))

ax1 = plt.subplot2grid((2,3),(0,0))
ax2 = plt.subplot2grid((2,3),(0,1))
ax3 = plt.subplot2grid((2,3),(0,2))
ax4 = plt.subplot2grid((2,3),(1,0), colspan=2)
ax5 = plt.subplot2grid((2,3),(1,2))

# ====================================================
# Animation Plot
# ====================================================

ax1.set_xlim(-2.5,2.5)
ax1.set_ylim(-1,1)

spring, = ax1.plot([],[],lw=3)
mass, = ax1.plot([],[],'ro',markersize=14)

ax1.set_title("Oscillator")

# ====================================================
# Phase Space
# ====================================================

ax2.plot(x,v,color="blue")
phase_point, = ax2.plot([],[],'ro')

ax2.set_title("Phase Space")
ax2.set_xlabel("Position")
ax2.set_ylabel("Velocity")

# ====================================================
# Energy
# ====================================================

ax3.plot(t,Energy,color="green")
energy_dot, = ax3.plot([],[],'ro')

ax3.set_title("Mechanical Energy")

# ====================================================
# FFT
# ====================================================

ax4.plot(xf[positive],yf[positive])

ax4.set_xlim(0,5)
ax4.set_title("Frequency Spectrum")

# ====================================================
# Poincare
# ====================================================

ax5.scatter(px,pv,s=12,color="black")
ax5.set_title("Poincare Section")

# ====================================================
# Animation
# ====================================================

def animate(i):

    xx = x[i]

    xs = np.linspace(0,xx,60)
    ys = 0.05*np.sin(25*np.pi*np.linspace(0,1,60))

    spring.set_data(xs,ys)
    mass.set_data([xx],[0])

    phase_point.set_data([x[i]],[v[i]])

    energy_dot.set_data([t[i]],[Energy[i]])

    return spring,mass,phase_point,energy_dot


ani = FuncAnimation(
    fig,
    animate,
    frames=N,
    interval=5,
    blit=True
)

plt.tight_layout()
plt.show()
