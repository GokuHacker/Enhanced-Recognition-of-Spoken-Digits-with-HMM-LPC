# Speech Digit Recognition System 🗣️🔢

This project implements a speech digit recognition system that recognizes spoken digits from 0 to 9. The system is built from scratch using Hidden Markov Models (HMM) and Linear Predictive Coding (LPC) for feature extraction. The project uses cepstral coefficients and the Durbin algorithm for HMM parameter estimation.

## Table of Contents 📚

- [Overview](#overview)
- [Features](#features)
- [Technologies Used](#technologies-used)
- [Project Structure](#project-structure)
- [Detailed Explanation](#detailed-explanation)
  - [Normalization and Feature Extraction](#normalization-and-feature-extraction)
  - [Universe Creation](#universe-creation)
  - [Codebook Generation](#codebook-generation)
  - [Hidden Markov Model Training](#hidden-markov-model-training)
  - [Recognition](#recognition)
- [Contributing](#contributing)
- [License](#license)

## Overview 🎯

The system processes speech data to recognize digits by:
- Normalizing the audio data.
- Extracting features using Linear Predictive Coding (LPC) and Cepstral Coefficients.
- Training a Hidden Markov Model (HMM) for each digit.
- Recognizing the spoken digit based on the trained models.

## Features ✨

- **Normalization:** The audio data is normalized to reduce the impact of noise and variations in recording conditions. 🎶
- **Feature Extraction:** LPC and Cepstral Coefficients are used to extract meaningful features from the speech signal. 📉
- **Codebook Generation:** The system uses the Linde-Buzo-Gray (LBG) algorithm to generate a codebook from the universe of training samples. 📚
- **HMM Training:** The system trains HMMs for each digit using the extracted features. 🤖
- **Digit Recognition:** The system recognizes the spoken digit by comparing the features of the input speech with the trained models. 🏷️

## Technologies Used 🛠️

- **C++:** Core programming language. 💻
- **HMM:** Hidden Markov Models for speech recognition. 🧠
- **LPC:** Linear Predictive Coding for feature extraction. 🔍
- **Windows API:** Used for file handling and system operations. 🪟

## Project Structure 🗂️

- **Digit_Recognition.cpp:** Main source code for the digit recognition system. 📜
- **234101043_universe.csv:** Universe file storing the cepstral coefficients of the training data. 📈
- **234101043_codebook.csv:** Codebook generated from the universe set. 📑

## Detailed Explanation 📖

### Normalization and Feature Extraction 🔧

- **Normalization:** The audio data is normalized by removing the DC offset and scaling the amplitude to a predefined range. 🎚️
- **Feature Extraction:** The system uses Linear Predictive Coding (LPC) to calculate the reflection coefficients, which are then transformed into Cepstral Coefficients. These coefficients represent the speech signal in a compact form. 📊

### Universe Creation 🌌

The universe is a collection of feature vectors extracted from training samples. The system processes multiple audio files, normalizes them, and extracts cepstral coefficients, which are stored in a CSV file. 🗃️

### Codebook Generation 📘

The Linde-Buzo-Gray (LBG) algorithm is applied to the universe to generate a codebook. The codebook is a quantized representation of the feature space, which reduces the complexity of HMM training. 🛠️

### Hidden Markov Model Training 📚

Each digit is modeled using a Hidden Markov Model (HMM). The system initializes the HMM parameters and uses the extracted features to train the models. The training process involves calculating forward and backward probabilities (Alpha and Beta) and updating the HMM parameters. 🔄

### Recognition 🏆

The recognition process involves calculating the probability of the observed sequence given each HMM model. The model with the highest probability is chosen as the recognized digit. 🎯

## Contributing

Feel free to open issues or submit pull requests to contribute to this project.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
