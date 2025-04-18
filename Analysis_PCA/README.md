# Hyperspectral Analysis - PCA

This repository contains tools and scripts for performing Principal Component Analysis (PCA) on hyperspectral data.

## Getting Started

### Prerequisites

- MATLAB with App Designer
- Hyperspectral data in an appropriate format (Check in the Supplementary)

### Instructions

1. Run **PCA_imaging.mlapp** in App Designer to open the initial window.
2. Press the button **Load** to upload the spectral hypercube.

### Usage

The application allows users to load hyperspectral data and perform PCA to reduce the dimensionality of the data for further analysis.

### Screenshots

Below are some labeled windows of the application to guide you:

#### Window 1
![Window1](https://github.com/MartinaRiva0300/HyperspectralAnalysis/blob/main/Analysis_PCA/supplementary_material_README/Window1.jpeg)
*If you want to add a description.*

#### Window 2
![Window2](https://github.com/MartinaRiva0300/HyperspectralAnalysis/blob/main/Analysis_PCA/supplementary_material_README/Window2.jpeg)

#### Window 2.1
![Window 2.1](https://github.com/MartinaRiva0300/HyperspectralAnalysis/blob/main/Analysis_PCA/supplementary_material_README/Window%202.1.jpeg)

#### Window 3
![Window3](https://github.com/MartinaRiva0300/HyperspectralAnalysis/blob/main/Analysis_PCA/supplementary_material_README/Window%203.jpeg)

## Contributing

If you wish to contribute to this project, please follow the standard procedures for forking the repository, making changes, and submitting a pull request.

## License (TO FILL)

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Special thanks to everyone who contributed to the development of this tool.

## Supplementary: Hyperspectral Data Format (Check if it is true)

To use the MATLAB scripts provided in this repository, your hyperspectral data should be in the following format:

1. **File Format**: The data should be stored in CSV format.
2. **Data Structure**:
    - Each row represents a different sample.
    - Each column represents the intensity value at a specific wavelength.
3. **Required Metadata**:
    - Sample ID: A unique identifier for each sample.
    - Acquisition Date: The date when the sample was acquired.
4. **Example Data Entry**:
    ```csv
    SampleID, AcquisitionDate, Wavelength1, Wavelength2, ..., WavelengthN
    S001, 2025-03-11, 0.123, 0.456, ..., 0.789
    S002, 2025-03-12, 0.321, 0.654, ..., 0.987
    ```
5. **Preprocessing Steps**:
    - Ensure all intensity values are normalized between 0 and 1.
    - Remove any rows with missing or corrupted data.

Please ensure your data follows this format before running the analysis scripts.
