# A Knowledge–Data Dual-Driven Framework for Predicting the Molecular Properties of Rechargeable Battery Electrolytes

## Introduction

This repository contains the code accompanying the paper "A Knowledge–Data Dual-Driven Framework for Predicting the Molecular Properties of Rechargeable Battery Electrolytes". The framework is designed to facilitate the prediction of molecular properties essential for the development of advanced rechargeable battery electrolytes. The code is organized into several modules to handle different aspects of the framework, including data handling, organization, interpretability, prediction, and application.

## Installation

To set up the environment for running the code, we provide a configuration file. You can create the conda environment by running the following command:

```bash
conda env create -f kpi_env.yml
```

This command will automatically install all the required packages and dependencies as specified in the `kpi_env.yml` file.

## Repository Structure

The repository is organized into the following directories:

### 1. Data

This directory contains the raw data used for the analysis. Ensure that your dataset is properly placed in this directory before running the code.

### 2. Organization

The `Organization` module is responsible for data organization and statistical analysis. It includes scripts to preprocess the raw data, perform necessary transformations, and conduct statistical evaluations.

### 3. Interpretability

The `Interpretability` module focuses on interpretability and knowledge discovery. It contains scripts for generating insights from the data, including feature importance analysis and SHAP value computation.

### 4. Prediction

The `Prediction` module is dedicated to knowledge-based molecular property prediction. It includes machine learning models and algorithms designed to predict various molecular properties based on the provided data.

### 5. Application

The `Application` module demonstrates the application of the framework in specific scenarios, such as molecular neighbor search and high-throughput virtual screening.

### 6. Utilities

The `Utilities` directory contains helper scripts and utility functions that support the main modules. These scripts are essential for performing common tasks and ensuring the smooth operation of the framework.

## Contributions

We welcome contributions from the community. If you encounter any issues or have suggestions for improvements, please open an issue or submit a pull request. We appreciate your feedback and collaboration.

## License

This project is licensed under the MIT License.

## Acknowledgements

We would like to thank the contributors and the research community for their valuable input and support in developing this framework.
