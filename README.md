# lyapunov-stability-fluid-dynamics
Tools for identifying Lyapunov exponents (LEs) and Covariant Lyapunov Vectors (CLVs) of chaotic flows (OpenFOAM + MATLAB)

# Lyapunov and Stability Analysis of Chaotic Flows

This repository contains the codebase developed for our study on the chaotic dynamics of flow past two side-by-side square cylinders. It provides tools to compute Lyapunov exponents (LEs), Covariant Lyapunov Vectors (CLVs) using OpenFOAM and MATLAB.

📝 Abstract

We investigate the dynamics of the flow past two side-by-side square cylinders at Reynolds number 200 and gap ratio 1. The flow is highly irregular due to the complex interaction between the flapping jet emanating from the gap and the vortices shed in the wake. We first perform Spectral Proper Orthogonal Decomposition (SPOD) to understand the flow characteristics. We then conduct Lyapunov stability analysis by linearizing the Navier-Stokes equations around the irregular base flow and find that it has two positive Lyapunov exponents. The Covariant Lyapunov Vectors (CLVs) are also computed; they characterise the unstable, neutral, and stable directions of the tangent space. Contours of the time-averaged CLVs reveal that the footprint of the leading CLV is in the near-wake, whereas the other CLVs peak further downstream, indicating distinct regions of instability.

SPOD of the two unstable CLVs is then employed to extract the dominant coherent structures and oscillation frequencies. For the leading CLV, the two dominant frequencies match closely with the prevalent frequencies of the spectrum of the drag coefficient, and correspond to instabilities due to vortex shedding and jet flapping. The second unstable CLV captures the subharmonic instability of the shedding frequency. Global linear stability analysis (LSA) of the time-averaged flow identifies a neutral eigenmode that resembles the leading SPOD mode of the first CLV, with very similar structure and frequency. However, while LSA predicts neutrality, Lyapunov analysis reveals that this direction is unstable, exposing further the inherent limitations of the global LSA when applied to chaotic flows.

🎯 Purpose

The purpose of this repository is to enable fluid dynamicists to perform Lyapunov stability analysis on a wide range of chaotic fluid flows.

📊 Analyses Included

- Lyapunov exponent (LE) computation of 2D chaotic flows
- Covariant Lyapunov Vector (CLV) extraction using Gram-Schmidt vectors
- SPOD of flow and of unstable CLVs

🧰 Software Requirements

- C++
- OpenFOAM v10
- MATLAB (tested on R2022b)


📚 How to Cite

If you use any part of this repository, please cite our paper:

(Preprint DOI / arXiv link coming soon)

📬 Contact

Feel free to reach out with questions or suggestions:

- 📧 s.sahu22@imperial.ac.uk

📜 License

This project is licensed under the MIT License – see the [LICENSE](./LICENSE) file for details.
