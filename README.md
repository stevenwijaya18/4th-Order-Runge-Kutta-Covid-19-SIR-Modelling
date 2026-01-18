# 4th-Order-Runge-Kutta-Covid-19-SIR-Modelling

This MATLAB project simulates the spread of COVID-19 using the SIR (Susceptible-Infected-Recovered) epidemiological model. It solves the system of Ordinary Differential Equations (ODEs) using the numerical 4th Order Runge-Kutta (RK4) method and compares the simulation results against observed data from a CSV file.

## 📋 Project Overview

The script performs the following operations:
1.  Data Ingestion: Reads real-world COVID-19 data (S, I, R populations) from a CSV file.
2.  Initialization: Sets initial conditions ($S_0$, $I_0$, $R_0$) and epidemiological parameters ($\beta$, $\gamma$).
3.  Numerical Solution: Implements the RK4 algorithm to solve the SIR differential equations over a specified time domain.
4.  Visualization: Plots the numerical model (continuous lines) against the observed data (discrete markers) for comparison.

## 🧮 Mathematical Model

The SIR model divides the population into three compartments:
* S: Susceptible
* I: Infected
* R: Recovered

The dynamics are defined by the following set of ODEs:

$$\frac{dS}{dt} = -\beta S I$$

$$\frac{dI}{dt} = \beta S I - \gamma I$$

$$\frac{dR}{dt} = \gamma I$$

Where:
* $\beta$ (beta): The infection rate.
* $\gamma$ (gamma): The recovery rate.

## 📂 File Requirements

To run this script, you must have the following file in your MATLAB directory:

### `Data COVID-19.csv`
The script uses `dlmread` with a semicolon delimiter. The file should adhere to this format:
* Delimiter: Semicolon (`;`)
* Header: 1 row (skipped by the script)
* Columns:
    1. Susceptible Data
    2. Infected Data
    3. Recovered Data

## ⚙️ Configuration & Parameters

You can tune the simulation by modifying the variables in the "Numerical Setup" section of the script:

| Parameter | Variable | Value in Code | Description |
| :--- | :--- | :--- | :--- |
| Initial Susceptible | `S0` | `1709` | Starting healthy population |
| Initial Infected | `I0` | `323` | Starting active cases |
| Initial Recovered | `R0` | `468` | Starting recovered cases |
| Infection Rate | `beta` | `0.000025` | Probability of transmitting disease |
| Recovery Rate | `gamma` | `0.05` | Rate at which infected recover |
| Time Span | `tn` | `150.0` | Total days to simulate |

## 🚀 Usage

1.  Ensure `Data COVID-19.csv` is in the same folder as your `.m` script.
2.  Open the script in MATLAB.
3.  Run the script (F5).
4.  Figure 1 will appear, displaying:
    * Markers (x): Real data points.
    * Lines (-): Model simulation.

## 📊 Plot Legend
* Red: Susceptible ($S$)
* Green: Infected ($I$)
* Blue: Recovered ($R$)