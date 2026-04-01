# C2DTD – CARBON-2D Topological Descriptor

This repository provides an implementation of the **CARBON-2D Topological Descriptor (C2DTD)**, a physically informed and interpretable representation for 2D carbon systems.

## 🧠 Overview

C2DTD combines three key structural components:

* **Local geometry** (coordination, bond lengths, angles)
* **Medium-range order** (radial distribution function – RDF)
* **Topology** (primitive ring statistics)

The descriptor is:

* Compact
* Physically interpretable
* Efficient for small datasets
* Suitable for machine learning applications

---

## ⚠️ Important Note on Interpretability

Although the descriptor is physically meaningful, the raw output (`descriptors.csv`) is **high-dimensional and not directly interpretable by inspection**.

To address this, this repository provides a dedicated script:

```bash
python interpretability.py
```

This script transforms the descriptor into **human-interpretable visual and statistical insights**, enabling:

* Understanding of dataset-level topology (e.g., dominance of hexagonal rings)
* Structure-by-structure analysis of ring distributions
* Detection of disorder, defects, and reconstruction patterns
* Visualization of how topology varies across the dataset

---

## ⚙️ Requirements

* Python 3.x
* numpy
* pandas
* matplotlib
* networkx
* ASE (Atomic Simulation Environment)

Install with:

```bash
pip install numpy pandas matplotlib networkx ase
```

---

## 📂 Usage

### 1. Generate descriptors

```bash
python descriptor_c2dtd.py
```

Output:

```
descriptors.csv
```

---

### 2. Interpret the dataset (ESSENTIAL STEP)

```bash
python interpretability.py
```

Outputs:

* 📊 Average ring distribution (global topology)
* 🧩 Per-structure ring distributions
* 🌡️ Heatmap of structural diversity
* 📄 Statistical summary of the dataset

This step is **essential** to extract physical insight from the descriptor.

---

## 📊 Descriptor Structure

The final descriptor includes:

* Structural statistics (mean, std, min, max)
* Radial distribution histogram
* Ring fractions (3–22)

---

## 📄 License

This project is licensed under the MIT License.

---

## 📚 Citation

If you use this code, please cite:

*C2DTD – CARBON-2D Topological Descriptor (2025)*

---

## 🚀 Future Work

* Integration with ML pipelines
* Extension to other 2D materials
* Large-scale screening

