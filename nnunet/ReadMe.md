# Segmentation Model using nnU-Net v1

## Overview
We developed a U-Net-based segmentation model to facilitate the end-to-end application of models that require segmentation. Our approach is built upon **nnU-Net v1**, an automated framework for configuring U-Net-based segmentation methods.

The **nnU-Net** framework is designed to automatically adapt to the characteristics of a given dataset, achieving state-of-the-art performance across various biomedical segmentation tasks. It has demonstrated its efficacy on 23 public datasets used in international segmentation competitions.

## Model Details
Our segmentation pipeline utilizes the complete **nnU-Net stack**, including:
- **2D U-Net**: Captures fine-grained details in individual slices.
- **3D U-Net**: Exploits spatial correlations across the z-axis.
- **3D U-Net Low Resolution**: Efficiently processes large volumes with reduced computational demand.
- **Ensemble Model**: Combines the predictions from the different nnU-Net configurations to enhance performance and robustness.

### Training Process
- **1000 epochs** with **5-fold cross-validation**.
- **Dataset-specific adaptations** performed by nnU-Net.
- **Automated optimization** of network depth, width, and resolution settings based on dataset characteristics.

## Model Weights
We have open-sourced the trained model weights, which can be accessed [here](#). However, we **strongly recommend** training your own models following the nnU-Net training procedures to achieve optimal results tailored to your specific dataset.

## Installation & Usage
To use this model or train your own using nnU-Net v1, follow these steps:

### Prerequisites
Ensure you have the following dependencies installed:
- Python 3.7+
- PyTorch
- nnU-Net v1
- Necessary medical imaging libraries (SimpleITK, NumPy, etc.)

### nnU-Net Installation

pip install nnunetv1

For full installation details, refer to the official nnU-Net v1 repository: [https://github.com/MIC-DKFZ/nnUNet](https://github.com/MIC-DKFZ/nnUNet)

### Training with nnU-Net
To train a model using nnU-Net:

nnUNet_train 2d nnUNetTrainerV2 TaskXXX_FoldX
nnUNet_train 3d nnUNetTrainerV2 TaskXXX_FoldX
nnUNet_train 3d_lowres nnUNetTrainerV2 TaskXXX_FoldX

After training, you can ensemble the models:

nnUNet_find_best_configuration -m 2d 3d 3d_lowres -t TaskXXX

### Inference with Pretrained Weights
To run inference using the open-sourced weights:

nnUNet_predict -i /path/to/images -o /path/to/output -t TaskXXX -m 2d 3d 3d_lowres -f all

## Citation
If you use this model or nnU-Net in your work, please cite:

Isensee, Fabian, et al. "nnU-Net: a self-adapting framework for U-Net-based medical image segmentation." Nature Methods, 2021

