from PIL import Image
import numpy as np

# Carica l'immagine
img = Image.open('F1_fus.png').convert('RGB')
data = np.array(img)

# Definisci i colori target (RGB)
rosso = [255, 0, 0]
giallo = [255, 215, 0]
#nero = [0, 0, 0]
verde=[0, 128,0]

# Crea maschere booleane per ogni colore
mask_rosso = np.all(data == rosso, axis=-1)
mask_giallo = np.all(data == giallo, axis=-1)
#mask_nero = np.all(data == nero, axis=-1)
mask_verde = np.all(data == verde, axis=-1)

# Conta i pixel per ciascun colore
count_rosso = np.sum(mask_rosso)
count_giallo = np.sum(mask_giallo)
#count_nero = np.sum(mask_nero)
count_verde = np.sum(mask_verde)


# Calcola i totali e le percentuali
total = count_rosso + count_giallo + count_verde
perc_rosso = count_rosso / total * 100
perc_giallo = count_giallo / total * 100
#perc_nero = count_nero / total * 100
perc_verde = count_verde / total * 100


# Stampa i risultati
print(f"Rosso: {perc_rosso:.2f}%")
print(f"Giallo: {perc_giallo:.2f}%")
#print(f"Nero: {perc_nero:.2f}%")
print(f"Verde: {perc_verde:.2f}%")
print(f"Total area: {total}")
