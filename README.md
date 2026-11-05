# RDKit Molecular Generator Overview

A Python tool for generating 2D visualisations of molecular libraries using RDKit.

The script takes a path to a CSV file containing compound identifiers and SMILES representations, generates molecular structures, and groups compounds based on structural similarity.

## 

## Features

* Import molecular structures from CSV files.
* Generate molecular structures from SMILES.
* Compounds clustered using Butina clustering algorithm using MACCS molecular fingerprints and Tanimoto similarity.
* Optionally align molecules to a reference SMILES.
* Generate 2D grids of molecular structures.
* Export molecular structures to SDF format.

## 

## Technologies

* Python 3.11
* RDKit
* NumPy
* Matplotlib

## 

## Input

A CSV file containing two columns:

* **ID**: Compound identifier.
* **SMILES**: Molecular SMILES representation.

## 

## Output

* 2D molecular structure visualisations.
* Optional SDF file containing clustered compounds.

## 

## Installation

Create the Conda environment:

conda env create -f mol\_gen\_env



Activate it:
conda activate mol\_gen\_env

