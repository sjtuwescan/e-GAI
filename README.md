# e-GAI: e-Value-based Generalized \u03B1-Investing for Online FDR Control

This repository implements the **e-GAI** and **p-RAI** methods described in our paper:

> **e-GAI: e-value-based Generalized \u03B1-Investing for Online False Discovery Rate Control**

## 📁 Repository Structure

```
├── Simulation/
│   ├── e-GAI-RAI.R        # Main driver for simulation experiments
│   └── mem-e-GAI.R        # Memory-efficient FDR control simulations
│
├── functions/
│   ├── eLORD.R            # Implements the e-LORD algorithm
│   ├── eSAFFRON.R         # Implements the e-SAFFRON algorithm
│   ├── pL-RAI.R           # Implements the pL-RAI procedure
│   ├── pS-RAI.R           # Implements the pS-RAI procedure
│   └── README.md          # Algorithm-specific notes & references
│
├── Real_data.Rmd          # Analysis of real data (Section 5.2 & 5.3)
└── LICENSE                # MIT License
```

## 🚀 Usage

### 1. Functions (Core Algorithms)

Each script in the `functions/` directory implements one of the online FDR control algorithms:

* **eLORD.R**: e-LORD
* **eSAFFRON.R**: e-SAFFRON
* **pL-RAI.R**: pL-RAI
* **pS-RAI.R**: pS-RAI

```

### 2. Simulation

All simulation experiments from Section 5 and appendix of the paper are orchestrated in `Simulation/e-GAI-RAI.R` and `Simulation/mem-e-GAI.R`:

```bash
Rscript Simulation/e-GAI-RAI.R       # Full simulation with default settings
Rscript Simulation/mem-e-GAI.R       # Memory FDR control mode
```

Simulation outputs (FDR, power) will be saved to `output/` by default.

### 3. Real Data Analysis 

The `Real_data.Rmd` notebook reproduces the analyses from Sections 5.2 and 5.3:

1. Open the notebook in RStudio or Jupyter:

   ```r
   rmarkdown::render("Real_data.Rmd")
   ```
2. Follow the narrative to load data, apply the algorithms, and visualize results.

## 📖 References

* **Paper**: [e-GAI: e-value-based Generalized $\alpha$-Investing for Online FDR Control](link-to-pdf)
* **Supplementary**: See algorithms described in `functions/` for technical details.

## 📝 License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

*Last updated: 2025-05-26*

