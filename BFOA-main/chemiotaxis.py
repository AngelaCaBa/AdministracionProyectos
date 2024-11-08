import math
import random
import numpy as np
from bacteria import bacteria

class chemiotaxis():
    def __init__(self, max_population=30):
        self.parcialNFE = 0
        self.max_population = max_population

    def compute_cell_interaction(self, bacteria, poblacion, d, w):
        total = 0.0
        for other in poblacion:
            diff = (bacteria.blosumScore - other.blosumScore) ** 2.0
            total += d * math.exp(w * diff)
        return total

    def attract_repel(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
        attract = self.compute_cell_interaction(bacteria, poblacion, -d_attr, -w_attr)
        repel = self.compute_cell_interaction(bacteria, poblacion, h_rep, -w_rep)
        return attract + repel

    def chemio(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
        if bacteria.NFE < 5:  # Inicialización de NFE en 5
            bacteria.NFE = 5

        previous_fitness = bacteria.fitness
        bacteria.interaction = self.attract_repel(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
        bacteria.fitness = bacteria.blosumScore + bacteria.interaction

        # Solo incrementar NFE si la mejora en fitness es significativa
        if bacteria.fitness > previous_fitness * 1.05:  # Incremento de NFE solo si mejora más del 5%
            bacteria.NFE += 0.2  # Incremento muy lento en caso de mejora
        else:
            bacteria.NFE += 0.05  # Incremento aún más lento sin mejora significativa

    def doChemioTaxis(self, poblacion, d_attr, w_attr, h_rep, w_rep):
        self.parcialNFE = 0
        for bacteria in poblacion:
            self.chemio(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
            self.parcialNFE += min(0.2, bacteria.NFE / 10)  # Control estricto de incremento en NFE

        # Mantener la población controlada
        if len(poblacion) > self.max_population:
            poblacion[:] = poblacion[:self.max_population]

    def eliminarClonar(self, path, poblacion):
        poblacion.sort(key=lambda x: x.fitness)
        eliminar = int(len(poblacion) * 0.6)
        for _ in range(eliminar):
            del poblacion[0]

        clones = self.clonacion(path, poblacion)
        for clone in clones:
            if len(poblacion) < self.max_population:
                poblacion.append(clone)

    def clonacion(self, path, poblacion):
        poblacionClones = []
        best = max(poblacion, key=lambda x: x.fitness if not np.isnan(x.fitness) else float('-inf'))

        for bacteria in poblacion:
            if np.isnan(bacteria.fitness):
                bacteria.fitness = 0

            # Mutación aún más baja
            mutacion = max(1, int((best.fitness - bacteria.fitness) / 25))

            newBacteria = bacteria.clonar(path)
            newBacteria.tumboNado(mutacion)
            newBacteria.autoEvalua()
            poblacionClones.append(newBacteria)

        return poblacionClones

    def randomBacteria(self, path):
        bact = bacteria(path)
        bact.tumboNado(random.randint(1, 1))  # Limitamos la variación de mutación
        return bact

    def insertRamdomBacterias(self, path, num, poblacion):
        for _ in range(num):
            if len(poblacion) < self.max_population:
                poblacion.append(self.randomBacteria(path))
            if len(poblacion) > self.max_population:
                poblacion.sort(key=lambda x: x.fitness)
                del poblacion[0]
