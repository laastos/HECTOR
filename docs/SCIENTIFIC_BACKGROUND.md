# Scientific Background

This document provides the scientific context and methodology underlying the HECTOR algorithm and its application to de novo protein binder design.

---

## Table of Contents

1. [Introduction](#introduction)
2. [The Binder Design Challenge](#the-binder-design-challenge)
3. [Complementarity-First Approach](#complementarity-first-approach)
4. [Surface Fingerprinting Principles](#surface-fingerprinting-principles)
5. [Experimental Validation](#experimental-validation)
6. [Comparison with Other Methods](#comparison-with-other-methods)
7. [References](#references)

---

## Introduction

Protein-protein interactions are central to virtually all biological processes. The ability to design proteins that bind specifically to predetermined targets has immense value for:

- **Therapeutics**: Blocking disease-related protein interactions
- **Diagnostics**: Detecting biomarkers with high specificity
- **Research tools**: Probing cellular mechanisms
- **Biotechnology**: Creating novel functional proteins

HECTOR (High-throughput Epitope Complementarity-based Target-Oriented Ranking) was developed to address the challenge of designing binders capable of targeting arbitrarily selected epitopes without relying on template structures from known binding partners.

### Key Publication

> **"A Complementarity-Based Approach to De Novo Binder Design"**
> Maksymenko K, Hatskovska V, Coles M, et al.
> *Advanced Science*, 2025, 12, e02015

---

## The Binder Design Challenge

### Requirements for Ideal Binders

An effective protein binder must simultaneously achieve:

1. **High binding affinity**: Strong interaction with the target
2. **High specificity**: Minimal off-target binding
3. **Stability**: Folded and functional in presence/absence of target

These requirements often conflict - optimizing one property may compromise others.

### Traditional Approaches

#### Template-Based Design
- Uses structural information from known binders
- Key interaction residues incorporated into compatible scaffolds
- Limited to targets with known binding partners

#### De Novo Backbone Design
- Builds backbones around natural binding motifs
- Covers narrow range of targetable surfaces
- Computationally expensive

#### Machine Learning Methods
- Deep neural networks for complementarity inference
- Diffusion models for backbone generation
- Confined to patterns learned from training data
- Success rates typically <1% without extensive screening

### HECTOR's Innovation

HECTOR introduces a **training-free** approach that:
- Relies only on geometric principles and force field parameters
- Generalizes to any target without learned patterns
- Achieves high success rates (>25% nanomolar binders)

---

## Complementarity-First Approach

### Rationale

Shape complementarity is a primary driver of protein binding:
- Maximizes van der Waals contacts
- Excludes solvent from interface
- Provides geometric specificity

HECTOR prioritizes shape complementarity as the first optimization criterion, followed by sequence design for favorable interactions.

### Two-Stage Design Pipeline

```
Stage 1: HECTOR Docking
├── Identify scaffolds with complementary surfaces
├── Filter by geometric criteria
└── Generate initial complex models

Stage 2: Interface Design
├── Optimize scaffold residues for target interactions
├── Rosetta/Damietta sequence design
├── MD-based validation
└── Select top candidates for experimental testing
```

### Advantages of Decoupling

1. **Reduced computational cost**: Geometry-based search is fast
2. **Broader scaffold diversity**: Not limited to common binder folds
3. **Training-free**: No dependence on available structure databases

---

## Surface Fingerprinting Principles

### Spin Images

HECTOR fingerprints are based on **spin images**, a well-established shape descriptor from computer vision.

#### Original Formulation (Johnson & Hebert, 1999)

For an oriented point O with position p and unit normal n, the spin map S_O projects 3D points into a 2D representation:

```
S_O: R^3 -> R^2
S_O(x) -> (alpha, beta)

alpha = sqrt(||x-p||^2 - (n.(x-p))^2)   [radial distance]
beta = n.(x-p)                           [axial distance]
```

#### HECTOR Adaptation

Key modifications for protein surface fingerprinting:

1. **Cylindrical basis projection**: Fixed-size 2D maps
2. **Angular integration**: Rotation invariance
3. **Invertibility**: Trivial complementarity transformation
4. **Sigmoid fading**: Reduced edge contributions

### Invertibility Property

A crucial insight enabling HECTOR's efficiency:

- A **forward map** describes a receptor surface
- An **inverse map** describes the ideal complementary ligand
- Inversion requires only a single transformation (flip along axial dimension)

```python
# Forward mapping (receptor): k = +1
# Inverse mapping (ligand): k = -1

inverse_map = forward_map[::-1, :]  # Simple flip
```

### R-Factor Calculation

The R-factor quantifies dissimilarity between fingerprints:

```
R = sum((2*F*I + C) / (F^2 + I^2 + C)) / N
```

- Lower R-factor indicates higher complementarity
- Cutoff of -0.82 provides good enrichment
- Sub-microsecond evaluation on GPU enables massive screening

### Database Search Strategy

Two-vs-all search with multiple query patches:

1. Select 2-3 patches at target epitope
2. Independently search for complementary subject patches
3. Filter by inter-patch distance constraints
4. Rank by combined R-factor and structural RMSD

---

## Experimental Validation

### VEGF Binders (Sam and Sima)

**Target**: Vascular Endothelial Growth Factor (VEGF) receptor-binding site

**Results**:
| Design | Scaffold | Kd (nM) | Tm (deg C) |
|--------|----------|---------|------------|
| Sam0.7 | 1OH0 (ketosteroid isomerase) | 190 | 63 |
| Sima3.2 | 1PM1 (nitrophorin) | 14 | 65 |

**Validation**:
- Crystal structures match design models (RMSD 1.7-1.8 A)
- Competes with bevacizumab for VEGF binding
- Reduces VEGF-dependent cell survival in vitro
- Tumor-inhibiting activity in zebrafish xenografts

### IL-7Ralpha Binders (des01-des08)

**Target**: Interleukin-7 receptor alpha at two epitopes

**Results**:
| Design | Scaffold | Target Site | Kd (nM) |
|--------|----------|-------------|---------|
| des01 | 5NLC (TIM barrel) | Site 1 | 26 |
| des03 | 6B8F (four-helix bundle) | Site 1 | 1.4 |
| des06 | 6B8F | Site 1 | - |
| des07 | 6YUD (Rossmann fold) | Site 2 | 20 |

**Validation**:
- 6/8 designs bind IL-7Ralpha
- 3 designs inhibit IL-7 signaling in reporter assay
- Designs from 3 different scaffold folds succeeded

### Success Rate Analysis

| Metric | Value |
|--------|-------|
| VEGF binders tested | 16 |
| VEGF nanomolar binders | 4 (25%) |
| IL-7Ralpha designs tested | 8 |
| IL-7Ralpha binders identified | 6 (75%) |
| IL-7Ralpha nanomolar binders | 6 (75%) |

---

## Comparison with Other Methods

### Template-Based Methods

| Method | Advantage | Limitation |
|--------|-----------|------------|
| Hotspot grafting | High affinity achievable | Limited targets |
| Natural motif transplant | Validated binding modes | Narrow epitope coverage |
| **HECTOR** | Any epitope targetable | Requires MD validation |

### Machine Learning Methods

| Method | Training Data | Success Rate | Screening Required |
|--------|---------------|--------------|-------------------|
| SurfaceID | PDB interfaces | <1% | Thousands |
| RFdiffusion | PDB structures | ~5-10% | Hundreds |
| **HECTOR** | None (training-free) | 25-75% | Tens |

### Key Differentiators

1. **No training data dependency**
   - HECTOR uses only force field parameters
   - Generalizes to novel chemotypes

2. **Scaffold diversity**
   - Not limited to common binder folds
   - Can utilize any PDB structure

3. **Explicit geometric optimization**
   - Prioritizes shape complementarity
   - Transparent scoring (R-factor)

4. **High success rate**
   - Fewer candidates needed for testing
   - Cost-effective experimental validation

---

## Biological Applications

### Therapeutic Targets

HECTOR is particularly suited for:

1. **Blocking protein-protein interactions**
   - Growth factor-receptor binding (VEGF/VEGFR)
   - Cytokine signaling (IL-7/IL-7R)
   - Immune checkpoints

2. **Targeting novel epitopes**
   - Sites without known binders
   - Epitopes inaccessible to antibodies

3. **Alternative to antibodies**
   - Single-domain proteins
   - No glycosylation required
   - E. coli expression compatible

### Design Considerations

For optimal results:

1. **Epitope selection**
   - Concave surfaces favor binding
   - Avoid highly flexible regions
   - Consider glycosylation sites

2. **Scaffold properties**
   - Small, stable domains preferred
   - Expression system compatibility
   - Minimal post-translational modifications

3. **Interface design**
   - Optimize polar contacts
   - Balance hydrophobic/hydrophilic residues
   - Validate with MD simulations

---

## Future Directions

### Algorithm Improvements

1. **Autoencoder compression**
   - 128-element vector fingerprints
   - Faster pre-filtering

2. **GPU-native implementation**
   - Full pipeline acceleration
   - Real-time database screening

3. **Multi-patch queries**
   - Beyond two-vs-all search
   - Improved specificity

### Applications

1. **AlphaFold database integration**
   - Millions of predicted structures
   - Expanded scaffold diversity

2. **Semi-humanized binders**
   - Human proteome scaffolds
   - Reduced immunogenicity

3. **Function transplantation**
   - Binding site grafting
   - Enzyme active site design

---

## References

### Primary Publication

1. Maksymenko K, Hatskovska V, Coles M, et al. A Complementarity-Based Approach to De Novo Binder Design. *Adv Sci*. 2025;12:e02015.

### Spin Image Methodology

2. Johnson AE, Hebert M. Using spin images for efficient object recognition in cluttered 3D scenes. *IEEE Trans Pattern Anal Mach Intell*. 1999;21:433-449.

### Structural Alignment

3. Kabsch W. A solution for the best rotation to relate two sets of vectors. *Acta Crystallogr A*. 1976;32:922-923.

### Related Binder Design Methods

4. Cao L, et al. De novo design of picomolar SARS-CoV-2 miniprotein inhibitors. *Science*. 2020;370:426-431.

5. Watson JL, et al. De novo design of protein structure and function with RFdiffusion. *Nature*. 2023;620:1089-1100.

6. Gainza P, et al. De novo design of protein interactions with learned surface fingerprints. *Nature*. 2023;617:176-184.

### Therapeutic Targets

7. Ferrara N, Adamis AP. Ten years of anti-vascular endothelial growth factor therapy. *Nat Rev Drug Discov*. 2016;15:385-403.

8. Markovic I, Savvides SN. Modulation of signaling mediated by TSLP and IL-7 in inflammation, autoimmune diseases, and cancer. *Front Immunol*. 2020;11:1557.

### Protein Design Software

9. Fleishman SJ, et al. RosettaScripts: a scripting language interface to the Rosetta macromolecular modeling suite. *PLoS One*. 2011;6:e20161.

10. Maksymenko K, et al. Damietta: A runtime-efficient tensorized protein design platform. *Cell Rep Methods*. 2023;3:100560.

---

## Data Availability

- **HECTOR software**: https://doi.org/10.24433/CO.9243108.v1
- **Crystal structures**: PDB 8BL5, 8BL9
- **MD trajectories**: https://doi.org/10.5281/zenodo.14028991
