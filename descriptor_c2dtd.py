import os
import re
import csv
import numpy as np
import networkx as nx

from ase.io import read
from ase.neighborlist import NeighborList, neighbor_list
from ase.build import make_supercell

# ================================================================
# ========================= PARAMETERS ===========================
# ================================================================

STRUCTURE_FOLDER = "../structures"
OUTPUT_CSV = "descriptors.csv"

# -------- STRUCTURAL DESCRIPTOR (PHYSICAL, NO FFT)
CUTOFF = 1.8
MAX_NEIGHBORS = 6
DENSITY_SIGMA = 2.0
RDF_BINS = 40
RDF_MAX = 5.0

# -------- RING DESCRIPTOR (CORRECTED)
RING_CUTOFF = 1.8
RING_RANGE = list(range(3, 23))

# ================================================================
# ======================= AUXILIARY FUNCTIONS ====================
# ================================================================

def ordenar_arquivos_cif(lista):
    def extrair_numero(nome):
        m = re.search(r"(\d+)", nome)
        return int(m.group(1)) if m else float("inf")
    return sorted(lista, key=lambda x: (extrair_numero(x), x))


# ================================================================
# ========== PHYSICAL STRUCTURAL DESCRIPTOR (NO FFT)
# ================================================================

def descritor_estrutural_fisico(cif_path):
    atoms = read(cif_path)
    atoms.pbc = [True, True, False]

    positions = atoms.get_positions()
    N = len(atoms)

    nl = NeighborList([CUTOFF/2]*N, self_interaction=False, bothways=True)
    nl.update(atoms)

    coord_numbers = []
    local_densities = []
    bond_lengths = []
    angles = []
    all_distances = []

    for i in range(N):
        idxs, offsets = nl.get_neighbors(i)

        if len(idxs) == 0:
            continue

        vecs = positions[idxs] + np.dot(offsets, atoms.cell) - positions[i]
        d = np.linalg.norm(vecs, axis=1)

        coord_numbers.append(len(d))
        bond_lengths.extend(d.tolist())
        all_distances.extend(d.tolist())

        # Gaussian local density
        density = np.sum(np.exp(-(d**2)/(2*DENSITY_SIGMA**2)))
        local_densities.append(density)

        # Local angles
        k = min(len(vecs), MAX_NEIGHBORS)
        v = vecs[:k]

        for a in range(len(v)):
            for b in range(a+1, len(v)):
                va = v[a] / np.linalg.norm(v[a])
                vb = v[b] / np.linalg.norm(v[b])
                ang = np.degrees(
                    np.arccos(np.clip(np.dot(va, vb), -1, 1))
                )
                angles.append(ang)

    # Robust statistics (invariant)
    def stats(x):
        if len(x) == 0:
            return [0,0,0,0]
        x = np.array(x)
        return [
            np.mean(x),
            np.std(x),
            np.min(x),
            np.max(x)
        ]

    coord_stats = stats(coord_numbers)
    density_stats = stats(local_densities)
    bond_stats = stats(bond_lengths)
    angle_stats = stats(angles)

    # Compact RDF (very important for energy)
    hist, _ = np.histogram(
        all_distances,
        bins=RDF_BINS,
        range=(0, RDF_MAX),
        density=True
    )

    return np.concatenate([
        coord_stats,
        density_stats,
        bond_stats,
        angle_stats,
        hist
    ])


# ================================================================
# ===== PRIMITIVE RING FUNCTIONS (NO FALSE CYCLES)
# ================================================================

def _canonical_cycle(cycle):
    if cycle[0] == cycle[-1]:
        cycle = cycle[:-1]
    n = len(cycle)
    rotations = [tuple(cycle[i:] + cycle[:i]) for i in range(n)]
    rev = list(reversed(cycle))
    rotations_rev = [tuple(rev[i:] + rev[:i]) for i in range(n)]
    return min(rotations + rotations_rev)


def _is_chordless(G, cycle):
    k = len(cycle)
    cycle_set = set(cycle)
    for i in range(k):
        u = cycle[i]
        for v in G.neighbors(u):
            if v in cycle_set:
                j = cycle.index(v)
                if abs(i - j) not in (1, k-1):
                    return False
    return True


def primitive_rings(G, max_size=10):
    rings = set()
    edges = list(G.edges())

    for u, v in edges:
        if not G.has_edge(u, v):
            continue

        G.remove_edge(u, v)
        try:
            length = nx.shortest_path_length(G, u, v)
            if length + 1 > max_size:
                G.add_edge(u, v)
                continue

            for path in nx.all_shortest_paths(G, u, v):
                cycle = path + [u]
                canon = _canonical_cycle(cycle)
                if 3 <= len(canon) <= max_size:
                    if _is_chordless(G, list(canon)):
                        rings.add(canon)

        except nx.NetworkXNoPath:
            pass

        G.add_edge(u, v)

    return [list(r) for r in rings]


# ================================================================
# ========== RING DESCRIPTOR (PHYSICAL AND CORRECT)
# ================================================================

def contagem_aneis_carbono_2D(cif_path):
    atoms = read(cif_path)
    atoms.pbc = [True, True, False]

    # Moderate supercell (avoids periodic artifacts)
    P = [[3,0,0],[0,3,0],[0,0,1]]
    sc = make_supercell(atoms, P)

    i, j, d = neighbor_list("ijd", sc, cutoff=RING_CUTOFF)

    G = nx.Graph()
    for a, b in zip(i, j):
        a = int(a)
        b = int(b)
        if a != b:
            G.add_edge(a, b)

    if G.number_of_nodes() < 3:
        return np.zeros(len(RING_RANGE))

    try:
        ciclos = primitive_rings(G, max_size=max(RING_RANGE))
    except Exception:
        return np.zeros(len(RING_RANGE))

    contagem = {n: 0 for n in RING_RANGE}

    for c in ciclos:
        n = len(c)
        if n in contagem:
            contagem[n] += 1

    total = sum(contagem.values())

    if total == 0:
        return np.zeros(len(RING_RANGE))

    return np.array([contagem[n] / total for n in RING_RANGE])


# ================================================================
# ============================= MAIN =============================
# ================================================================

def main():

    arquivos = ordenar_arquivos_cif(
        [f for f in os.listdir(STRUCTURE_FOLDER) if f.endswith(".cif")]
    )

    linhas = []

    for f in arquivos:
        path = os.path.join(STRUCTURE_FOLDER, f)
        print("Processing:", f)

        desc_struct = descritor_estrutural_fisico(path)
        desc_rings = contagem_aneis_carbono_2D(path)

        desc_final = np.concatenate([desc_struct, desc_rings])
        linhas.append([f] + desc_final.tolist())

    with open(OUTPUT_CSV, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)

        header = (
            ["structure"]
            + [f"struct_{i}" for i in range(len(desc_struct))]
            + [f"ring_{n}" for n in RING_RANGE]
        )

        writer.writerow(header)
        writer.writerows(linhas)

    print("\n✅ PHYSICAL DESCRIPTOR (NO FFT) + RINGS GENERATED!")
    print("File:", OUTPUT_CSV)
    print("Descriptor dimension:", len(linhas[0]) - 1)


if __name__ == "__main__":
    main()
