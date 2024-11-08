from fastaReader import fastaReader
import random
import numpy as np
import copy
from evaluadorBlosum import evaluadorBlosum

class bacteria:
    def __init__(self, path):
        self.matrix = fastaReader(path)
        self.blosum = evaluadorBlosum()  # Instancia del evaluador BLOSUM
        self.blosumScore = 0
        self.fitness = 0
        self.interaction = 0
        self.NFE = 0
        self.autoEvalua()  # Evalúa automáticamente al crear una bacteria

    def autoEvalua(self):
        """Calcula el blosumScore y el fitness inicial de la bacteria."""
        self.blosumScore = self.calcula_blosum_score()
        self.fitness = self.blosumScore + self.interaction  # Incluye interacción en el cálculo
        self.NFE += 1  # Incrementa el NFE por esta evaluación

    def calcula_blosum_score(self):
        """Calcula el score de BLOSUM para la bacteria actual basado en las secuencias."""
        score = 0
        for seqA in self.matrix.seqs:
            for seqB in self.matrix.seqs:
                if seqA != seqB:
                    score += sum(self.blosum.getScore(a, b) for a, b in zip(seqA, seqB))
        return score

    def showGenome(self):
        """Muestra el genoma (secuencias) de la bacteria."""
        for seq in self.matrix.seqs:
            print(seq)

    def clonar(self, path):
        """Clona la bacteria actual, copiando las secuencias y calculando el fitness del clon."""
        newBacteria = bacteria(path)
        newBacteria.matrix.seqs = np.array(copy.deepcopy(self.matrix.seqs))
        newBacteria.autoEvalua()  # Calcula el fitness del clon
        return newBacteria

    def tumboNado(self, numGaps):
        """Inserta aleatoriamente un número de gaps en las secuencias de la bacteria."""
        self.cuadra()  # Ajusta las secuencias a la misma longitud
        matrixCopy = copy.deepcopy(self.matrix.seqs).tolist()  # Convierte a lista para modificar

        for _ in range(random.randint(0, numGaps)):  # Inserta un número aleatorio de gaps
            seqnum = random.randint(0, len(matrixCopy) - 1)  # Selecciona secuencia
            pos = random.randint(0, len(matrixCopy[0]))
            part1 = matrixCopy[seqnum][:pos]
            part2 = matrixCopy[seqnum][pos:]
            matrixCopy[seqnum] = "-".join([part1, part2])  # Inserta gap

        self.matrix.seqs = np.array(matrixCopy)  # Convierte de vuelta a numpy array
        self.cuadra()  # Ajusta nuevamente el tamaño
        self.limpiaColumnas()  # Limpia columnas de gaps innecesarios
        self.autoEvalua()  # Recalcula el fitness después de modificar la bacteria

    def cuadra(self):
        """Ajusta las secuencias de la bacteria para que tengan la misma longitud."""
        max_len = max(len(seq) for seq in self.matrix.seqs)
        for i in range(len(self.matrix.seqs)):
            self.matrix.seqs[i] = self.matrix.seqs[i].ljust(max_len, '-')  # Rellena con gaps

    def limpiaColumnas(self):
        """Elimina columnas con solo gaps."""
        columns = np.array([list(seq) for seq in self.matrix.seqs]).T  # Transpone las secuencias
        cleaned_columns = [col for col in columns if '-' in col and len(set(col)) > 1]
        self.matrix.seqs = np.array([''.join(col) for col in np.array(cleaned_columns).T])  # Reconstruye las secuencias
