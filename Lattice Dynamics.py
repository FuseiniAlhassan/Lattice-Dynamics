import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Output directory
FIG_DIR = "figures_lattice"
os.makedirs(FIG_DIR, exist_ok=True)

# Parameters
N = 20          # number of atoms
m = 1.0         # mass of each atom
k = 1.0         # spring constant
a = 1.0         # lattice spacing
dt = 0.05       # time step
n_steps = 400   # total time steps

# Initial positions and velocities

x0 = np.arange(N) * a
u = np.zeros(N)          # displacement
v = np.zeros(N)          # velocity

# Small initial perturbation for first atom
u[0] = 0.5

# Precompute dynamical matrix for normal modes

D = np.zeros((N,N))
for i in range(N):
    D[i,i] = 2*k/m
    D[i,(i-1)%N] = -k/m
    D[i,(i+1)%N] = -k/m

eigvals, eigvecs = np.linalg.eigh(D)
frequencies = np.sqrt(eigvals)

# Time evolution (Verlet integration)

displacements = []

for step in range(n_steps):
    # Compute forces
    F = np.zeros(N)
    for i in range(N):
        F[i] = -k*(2*u[i] - u[(i-1)%N] - u[(i+1)%N])
    # Update velocities and positions (velocity Verlet)
    v += 0.5*F*dt
    u += v*dt
    # Recompute forces
    F = np.zeros(N)
    for i in range(N):
        F[i] = -k*(2*u[i] - u[(i-1)%N] - u[(i+1)%N])
    v += 0.5*F*dt

    displacements.append(u.copy())

displacements = np.array(displacements)

# Plot atomic displacements at selected times

fig, ax = plt.subplots(figsize=(10,4))
for t in [0, 50, 100, 200, 300, 399]:
    ax.plot(x0, displacements[t], label=f"t={t*dt:.2f}")
ax.set_xlabel("Atom index")
ax.set_ylabel("Displacement")
ax.set_title("Atomic displacements over time")
ax.grid(True)
ax.legend()
fname = os.path.join(FIG_DIR, "displacements_over_time.png")
plt.savefig(fname, dpi=200)
plt.close(fig)
print("Saved figure:", fname)

# Phonon dispersion (normal modes)

fig, ax = plt.subplots(figsize=(8,4))
q = 2*np.pi*np.arange(N)/N
omega = frequencies
ax.plot(q, omega, 'o-')
ax.set_xlabel("Wavevector q")
ax.set_ylabel("Frequency ω")
ax.set_title("Phonon dispersion (1D chain)")
ax.grid(True)
fname = os.path.join(FIG_DIR, "phonon_dispersion.png")
plt.savefig(fname, dpi=200)
plt.close(fig)
print("Saved figure:", fname)

# Animation of lattice vibrations
fig, ax = plt.subplots(figsize=(10,4))
line, = ax.plot(x0, displacements[0], 'o-', lw=2)
ax.set_ylim(-1.2*np.max(np.abs(displacements)), 1.2*np.max(np.abs(displacements)))
ax.set_xlabel("Atom index")
ax.set_ylabel("Displacement")
ax.set_title("Lattice vibration animation")
ax.grid(True)

def animate(i):
    line.set_ydata(displacements[i])
    ax.set_title(f"Lattice vibration t={i*dt:.2f}")
    return line,

anim = FuncAnimation(fig, animate, frames=n_steps, interval=50, blit=True)
anim_fname = os.path.join(FIG_DIR, "lattice_vibration.gif")
anim.save(anim_fname, writer='pillow', fps=25)
plt.close(fig)
print("Saved animation GIF:", anim_fname)
