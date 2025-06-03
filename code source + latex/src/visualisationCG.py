import numpy as np
import matplotlib.pyplot as plt

# Charger les données depuis le fichier binaire
def load_phi(filename):
    with open(filename, 'rb') as f:
        # Lire les dimensions de la matrice
        nx, ny = np.fromfile(f, dtype=np.int32, count=2)
        print(f"Dimensions lues: nx = {nx}, ny = {ny}")
        
        # Lire les données de la matrice phi
        phi = np.fromfile(f, dtype=np.float64)
        print(f"Taille des données lues: {phi.size}")
        
        # Vérifier si la taille des données lues est correcte
        expected_size = nx * ny
        print(f"Taille attendue : {expected_size}")
        
        if phi.size != expected_size:
            print(f"Attention : La taille des données lues ({phi.size}) ne correspond pas aux dimensions attendues ({expected_size})")
            # Ajuster la taille des données pour correspondre à la taille attendue
            phi = phi[:expected_size]  # Utiliser seulement la taille attendue
        
        # Redimensionner les données avec les dimensions lues
        phi = phi.reshape((ny, nx))
    
    return phi

# Charger les données
phi = load_phi('solution_cg.dat')

# Visualisation avec contourf (niveau de contours de 50)
plt.figure(figsize=(6, 5))
plt.contourf(phi, levels=50, cmap="inferno_r")  # Utiliser "inferno" pour les couleurs
plt.colorbar(label="Valeur de Phi")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Solution de l'équation Phi - ΔPhi = S en 2D")
plt.show()
