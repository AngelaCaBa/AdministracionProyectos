from bacteria import bacteria
from chemiotaxis import chemiotaxis
import numpy

# Definir parámetros
poblacion = []
path = "C:/Users/angel/OneDrive - Universidad Autonoma de Coahuila/7 Semestre/Administracion de proyectos/BFOA-main/multiFasta.fasta"
numeroDeBacterias = 5
numRandomBacteria = 1
iteraciones = 10
tumbo = 1  # Número de gaps a insertar
nado = 3
chemio = chemiotaxis()

# Bacterias
veryBest = bacteria(path)  # Mejor bacteria
tempBacteria = bacteria(path)  # Bacteria temporal para validaciones
original = bacteria(path)  # Bacteria original sin gaps

# Global NFE
globalNFE = 0  # Número de evaluaciones de la función objetivo

# Parámetros de quimiotaxis
dAttr = 0.1  # 0.1
wAttr = 0.2  # 0.2
hRep = dAttr
wRep = 10  # 10

def clonaBest(veryBest, best):
    """Clonar la mejor bacteria a 'veryBest'."""
    if best.fitness > veryBest.fitness:
        veryBest.matrix.seqs = numpy.array(best.matrix.seqs)
        veryBest.blosumScore = best.blosumScore
        veryBest.fitness = best.fitness
        veryBest.interaction = best.interaction

def validaSecuencias(path, veryBest):
    """Validar que las secuencias de 'veryBest' coincidan con las originales."""
    # Clonar a veryBest en tempBacteria
    tempBacteria.matrix.seqs = numpy.array(veryBest.matrix.seqs)
    
    # Eliminar gaps de las secuencias de tempBacteria
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-", "")
    
    # Validar que las secuencias originales coincidan con las secuencias de tempBacteria
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return

def optimizarQuimiotaxis(chemio, poblacion, dAttr, wAttr, hRep, wRep):
    """Optimización del proceso de quimiotaxis para usar menos NFE."""
    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)
    return chemio.parcialNFE

def eliminarClonar(chemio, path, poblacion):
    """Optimizar eliminación y clonación de bacterias para reducir NFE."""
    poblacion.sort(key=lambda x: x.fitness)
    
    # Eliminar el 20% con menor fitness
    for i in range(int(len(poblacion) * 0.2)):
        del poblacion[0]

    clones = chemio.clonacion(path, poblacion)
    for clone in clones:
        poblacion.append(clone)

def insertRamdomBacterias(chemio, path, numRandomBacteria, poblacion):
    """Insertar bacterias aleatorias optimizando la calidad de la población."""
    chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)

# Inicializar la población
for i in range(numeroDeBacterias):
    poblacion.append(bacteria(path))

# Ciclo de iteraciones
for _ in range(iteraciones):
    for bacteria in poblacion:
        bacteria.tumboNado(tumbo)
        bacteria.autoEvalua()  # Evaluar la bacteria
    globalNFE += optimizarQuimiotaxis(chemio, poblacion, dAttr, wAttr, hRep, wRep)  # Optimizar quimiotaxis

    # Actualizar la mejor bacteria encontrada
    best = max(poblacion, key=lambda x: x.fitness)
    if (veryBest == None) or (best.fitness > veryBest.fitness):
        clonaBest(veryBest, best)
    
    print(f"Interacción: {veryBest.interaction} | Fitness: {veryBest.fitness} | NFE: {globalNFE}")
    
    # Eliminar bacterias con bajo rendimiento y clonar
    eliminarClonar(chemio, path, poblacion)
    
    # Insertar bacterias aleatorias
    insertRamdomBacterias(chemio, path, numRandomBacteria, poblacion)

    print(f"Población: {len(poblacion)}")

# Mostrar el genoma de la mejor bacteria
veryBest.showGenome()

# Validar que las secuencias son correctas
validaSecuencias(path, veryBest)
