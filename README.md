# Radar Declutter

This repo is a Python package to correct radargrams captured from helicopter‐based echo soundings. It removes terrain‐induced noise by compensating for surface topography using a clutter modeling approach based on viewsheds (using [ArcGIS](https://pro.arcgis.com/)) and digital elevation data.

## Installation

The declutter code runs with Python 3. To run the code, the [ArcGIS Python API](https://developers.arcgis.com/python/latest/) is required. It is recommended to [create a conda environment](https://developers.arcgis.com/python/latest/guide/understanding-conda/#:~:text=providing%20more%20details.-,Conda%20Environments,versions%20of%20software%2C%20including%20Python) within your ArcGIS installation and run the model from there. You will also need `jupyter notebook` or `jupyter lab` to run the [example notebook](https://github.com/bearecinos/radar-declutter/blob/master/docs/Examples.ipynb).

## Quick Start / Usage

A simple demonstration is provided in the [example notebook](https://github.com/bearecinos/radar-declutter/blob/master/docs/Examples.ipynb), which includes:

- Processing terrain data
- Setting radargram parameters
- Generating synthetic radargrams from flight paths

See the [Wiki](https://github.com/bearecinos/radar-declutter/wiki) for the full workflow and module documentation.

## Citation

[![DOI](https://zenodo.org/badge/741457808.svg)](https://doi.org/10.5281/zenodo.15488953)