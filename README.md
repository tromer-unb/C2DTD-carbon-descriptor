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

### 2. Generate plots

```bash
python your_plot_script.py
```

Outputs:

* Average ring distribution
* Per-structure plots
* Heatmap
* Statistical summary

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
