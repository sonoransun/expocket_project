# expocket — Pre-miRNA DICER Cleavage Analysis & 3D Visualization

Analysis of how terminal nucleotide identity and structural features (YCR motif, overhang) control DICER cleavage site selection and accuracy on pre-miRNA substrates, with an interactive 3D visualization application.

## Project Overview

This repository contains two components:

1. **Analysis scripts** — Python scripts for processing high-throughput dicing assay data from a randomized pre-mir-324 library
2. **3D Viewer** — A GPU-accelerated desktop application for interactive 3D visualization of the sequence and cleavage data

```mermaid
graph TD
    A[Pre-mir-324 Library<br/>256 variants] --> B[High-Throughput<br/>Dicing Assay]
    B --> C[Raw Sequencing Counts]
    C --> D[Analysis Scripts]
    D --> E[Cleavage Accuracy<br/>& Efficiency Metrics]
    E --> F[2D Plots<br/>matplotlib / seaborn]
    E --> G[3D Viewer App<br/>pygfx / WebGPU]
```

## Analysis Scripts

| File | Description |
|------|-------------|
| `hsa_324_data_khoa.py` | Human DICER analysis of pre-mir-324 library (PNK/non-PNK conditions). Heatmaps, boxplots, reproducibility, DC21/DC22 accuracy, logo plots. |
| `dme_324lib_3repeat_khoa.py` | Fly Dicer-1 (DCR-1) analysis of pre-mir-324 library. Cross-species comparison with human DICER. |
| `YCR_end.py` | Analysis of YCR motif and 5'-nucleotide coordination in human pre-miRNAs (miRGeneDB data). |

### Experimental Design

```mermaid
graph LR
    subgraph "256 Variants"
        direction TB
        A["Group A (5'=A)<br/>64 variants"]
        T["Group T (5'=U)<br/>64 variants"]
        G["Group G (5'=G)<br/>64 variants"]
        C["Group C (5'=C)<br/>64 variants"]
    end

    subgraph "Per Variant"
        direction TB
        S["3' randomized<br/>3-nt sequence"] --> F["Fold &<br/>Structure"]
        F --> CL["DICER Cleavage"]
        CL --> DC20["DC20"]
        CL --> DC21["DC21"]
        CL --> DC22["DC22"]
        CL --> DC23["DC23"]
    end

    A --> S
    T --> S
    G --> S
    C --> S
```

### Key Metrics

- **Cleavage accuracy** — proportion of reads at a specific cleavage site relative to total reads for that variant
- **Positional efficiency** — log2-ratio of cleavage product RPM to control RPM at a specific position
- **Global efficiency** — log2-ratio of total cleavage product RPM to control RPM for a variant

## 3D Viewer Application

A standalone, cross-platform desktop application with two linked 3D panels for exploring the cleavage data interactively.

### Architecture

```mermaid
graph TB
    subgraph "Desktop Application"
        direction TB
        MW["Main Window<br/>(PySide6 / Qt6)"]

        subgraph "Left Panel"
            RNA["RNA Hairpin View<br/>3D molecular structure"]
        end

        subgraph "Right Panel"
            LS["Data Landscape<br/>3D scatter of 256 variants"]
        end

        subgraph "Controls"
            SB["Sidebar<br/>Cleavage site · Enzyme · Group filter"]
            IP["Info Panel<br/>Selected variant details"]
        end

        MW --> RNA
        MW --> LS
        MW --> SB
        MW --> IP
    end

    IC["Interaction Controller"] <--> RNA
    IC <--> LS
    SB --> IC
    IC --> IP

    subgraph "Rendering"
        PG["pygfx Scene Graph"] --> WG["wgpu-py<br/>(WebGPU)"]
        WG --> Metal["Metal<br/>(macOS)"]
        WG --> Vulkan["Vulkan<br/>(Linux)"]
        WG --> DX12["DX12<br/>(Windows)"]
    end

    RNA --> PG
    LS --> PG
```

### Panels & Interactions

```mermaid
sequenceDiagram
    participant User
    participant Landscape as Data Landscape
    participant Controller as Interaction Controller
    participant RNA as RNA Hairpin View
    participant Info as Info Panel

    User->>Landscape: Click variant sphere
    Landscape->>Controller: variant_clicked("G_ATC")
    Controller->>RNA: set_variant(G_ATC)
    Note over RNA: Rebuild hairpin<br/>with G_ATC structure
    Controller->>Landscape: highlight_variant("G_ATC")
    Note over Landscape: Glow on selected sphere
    Controller->>Info: update_variant("G_ATC")
    Note over Info: Show sequence, structure,<br/>cleavage data

    User->>Controller: Change to DC22
    Controller->>Landscape: Rebuild scatter (z = DC22 accuracy)
    Controller->>RNA: Update cleavage disc sizes
```

### Features

| Feature | Description |
|---------|-------------|
| **RNA 3D hairpin** | Nucleotide bases as colored spheres, backbone as spline tube, hydrogen bonds between paired bases |
| **Cleavage site markers** | Semi-transparent discs at DC20–DC23 positions, sized by cleavage accuracy |
| **Data landscape** | 256 variant spheres in a 4-group × 8×8 grid, height = accuracy, color = nucleotide group |
| **Cross-panel linking** | Selecting a variant in the landscape updates the RNA structure and info panel |
| **Enzyme toggle** | Switch between Human DICER and Fly DCR-1 datasets |
| **Group filtering** | Show/hide variants by 5' nucleotide group (A, T, G, C) |
| **Mock data** | Automatically generates realistic demo data when experimental files are unavailable |

### Running the Viewer

```bash
# Create and activate virtual environment
python3 -m venv .venv
source .venv/bin/activate    # macOS/Linux
# .venv\Scripts\activate     # Windows

# Install
pip install -e .

# Launch
python -m viewer
```

### Module Structure

```
viewer/
  __main__.py              Entry point
  app.py                   Main window (dual-panel layout)
  config.py                Colors, helix parameters, paths
  data/
    schema.py              VariantInfo, CleavageRecord, VariantDataset
    loader.py              Load TSV data files
    mock_data.py           Generate demo data
  rna3d/
    layout.py              Dot-bracket → 3D coordinates (A-form helix)
    scene.py               pygfx scene graph for RNA hairpin
    widgets.py             Qt widget with embedded pygfx canvas
  landscape/
    scene.py               pygfx scene graph for variant scatter
    widgets.py             Qt widget with embedded pygfx canvas
  interaction/
    controller.py          Cross-panel selection synchronization
  ui/
    sidebar.py             Filter & display controls
    info_panel.py          Variant detail display
```

## Dependencies

- **Python** ≥ 3.10
- **pygfx** + **wgpu** — GPU-accelerated 3D rendering via WebGPU
- **PySide6** — Qt6 desktop windowing
- **pandas**, **numpy**, **scipy** — data handling
- **matplotlib**, **seaborn** — 2D plotting (analysis scripts)
