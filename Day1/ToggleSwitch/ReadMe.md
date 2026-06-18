# Genetic Toggle Switch

## Biological Background

The **Genetic Toggle Switch** is a synthetic gene regulatory circuit engineered by Gardner et al. (*Nature*, 2000). It consists of two genes whose protein products mutually repress each other's transcription.

In the cell, **Gene 1** produces **Protein A**, which represses the transcription of Gene 2. At the same time, **Gene 2** produces **Protein B**, which represses the transcription of Gene 1.

```
Gene 1 ---[A]--->|  Gene 2
Gene 2 ---[B]--->|  Gene 1
```

This mutual repression creates a **bistable switch**: the system locks into one of two stable states:
- **State 1:** A high, B low (Gene 1 wins)
- **State 2:** A low, B high (Gene 2 wins)

Which state is reached depends on the initial conditions and stochastic fluctuations, making this an ideal model to study **bistability** and **noise-driven switching** with epimod.

---

## Mathematical Model

The dynamics are described by two ODEs using Hill-function repression:

$$\frac{dA}{dt} = \frac{\alpha_1}{1 + B^n} - \delta \cdot A$$

$$\frac{dB}{dt} = \frac{\alpha_2}{1 + A^m} - \delta \cdot B$$

where:

| Symbol | Meaning | Typical value |
|--------|---------|---------------|
| $\alpha_1$ | Max transcription rate of Gene 1 | 10–50 |
| $\alpha_2$ | Max transcription rate of Gene 2 | 10–50 |
| $n$ | Hill coefficient for B repressing Gene 1 | 2 |
| $m$ | Hill coefficient for A repressing Gene 2 | 2 |
| $\delta$ | Degradation/dilution rate (same for A and B) | 1.0 |

The bistable region exists when $n, m > 1$ and the production rates are roughly symmetric ($\alpha_1 \approx \alpha_2$).

---

## Steady States and Bistability

For symmetric parameters ($\alpha_1 = \alpha_2 = \alpha$, $n = m$), the two stable steady states are approximately:

$$A^* \approx \frac{\alpha}{\delta}, \quad B^* \approx 0 \qquad \text{(State 1)}$$
$$A^* \approx 0, \quad B^* \approx \frac{\alpha}{\delta} \qquad \text{(State 2)}$$

An unstable symmetric equilibrium $A^* = B^*$ separates the two basins of attraction.

---

## Running with epimod


**Suggested parameters to vary in `parameters_sensitivity.csv`:**

| Parameter | Distribution | Range |
|-----------|-------------|-------|
| $\alpha_1$ | uniform | [5, 50] |
| $\alpha_2$ | uniform | [5, 50] |
| $n$ | uniform | [1.5, 5] |
| $m$ | uniform | [1.5, 5] |

---

## Expected Results

### Deterministic trajectories
Starting from A=40, B=5 → system converges to **State 1** (A high).  
Starting from A=5, B=40 → system converges to **State 2** (B high).

### Stochastic simulations
With symmetric initial conditions and small populations (~10–20 molecules), individual runs switch randomly between states over time.

---

## References

Gardner T.S., Cantor C.R., Collins J.J. (2000).  
*Construction of a genetic toggle switch in Escherichia coli.*  
**Nature**, 403, 339–342. https://doi.org/10.1038/35002131
