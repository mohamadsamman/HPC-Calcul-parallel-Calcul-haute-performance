import numpy as np
import matplotlib.pyplot as plt

# Chargement des données
data = np.loadtxt("solution.dat")
x = data[:, 0]
y = data[:, 1]
phi = data[:, 2]

# Reshape en grille 2D
N = int(np.sqrt(len(x)))
X = x.reshape(N, N)
Y = y.reshape(N, N)
Z = phi.reshape(N, N)

# Affichage du champ
plt.figure(figsize=(6,5))
plt.contourf(X, Y, Z, levels=50, cmap="inferno")
plt.colorbar(label="Valeur de Phi")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Solution de l'équation Phi - ΔPhi = S en 2D")
plt.show()

