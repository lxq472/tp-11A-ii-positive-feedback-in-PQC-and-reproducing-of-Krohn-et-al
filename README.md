# Week 2 – Positive Feedback in PQC & Reproducing Krohn et al. JCI 2011

**Course:** Physics of Molecular Diseases, Prof. Ala Trusina, Niels Bohr Institute, 2020  
**Topic:** Protein Folding Diseases II — Protein Quality Control (PQC)  
**Question 2:** *"Add the effect of aggregates on PQC (Agg --|PQC--|Aggs: positive feedback). Can you reproduce observations by Krohn et al. JCI 2011 (delayed spike)?"*

---

## The biological story

In the previous simulation (Week 2, Part 1) we modelled PQC as a **negative feedback** that fights aggregates:

```
Aggregates  →  upregulate PQC  --|  Aggregates
```

But the cell cannot keep up indefinitely. **Insoluble aggregates** (plaques, fibrils) physically sequester and inhibit PQC components — chaperones, proteasomes, ABC transporters. This adds a **positive feedback loop** on top:

```
Insoluble Ai  --|  PQC  --|  monomer clearance
                            ↓
                        M rises
                            ↓
                  nucleation ∝ M^6 explodes
                            ↓
                        more As  →  more Ai   (loop closed)
```

This positive feedback is the central finding of **Krohn et al., JCI 2011**: the transporter ABCC1 is a key clearance component, and its loss (`ABCC1−/−` mice) **doubles the insoluble Aβ burden** and shifts the aggregation curve in a way consistent with a runaway positive feedback. Crucially, a mere **~11% increase in clearance capacity** (simulating `ABCB1−/−` vs `ABCC1−/−` transport kinetics) delays the onset of aggregation significantly — a highly non-linear response that is the hallmark of a system operating near a tipping point.

---

## What the model reproduces

The simulation quantitatively reproduces all three key observations from Krohn et al. fig. 4:

| Observation | Krohn fig. 4 | This model |
|-------------|--------------|------------|
| Soluble Aβ shows **delayed spike then decline** | fig. 4A/C | ✓ As peaks at t ≈ 19 wk then falls |
| Insoluble Aβ shows **delayed slow monotone growth** | fig. 4B/D | ✓ Ai grows continuously after t ≈ 10 wk |
| +11% clearance → **significant rightward shift** of spike | fig. 4C vs 4A | ✓ +11% cT → +3.9 wk delay (+21% relative) |
| No positive feedback → **no spike** (just monotone growth) | implicit | ✓ K→∞ removes spike entirely |

---

## Model equations

### Module 1 — monomer pool (Krohn)

$$\frac{dM}{dt} = p \;-\; \underbrace{\frac{c_T}{1+(A_i/K)^h}}_{\text{PQC, inhibited by }A_i} \cdot M \;-\; n_n \cdot k_n \cdot M^{n_n} \;-\; k_g \cdot M \cdot (A_s + A_i)$$

The second term is the **positive feedback node**: as $A_i$ grows, the Hill function reduces PQC efficiency, $M$ rises, and nucleation ($\propto M^{n_n}$) accelerates super-linearly.

### Module 2 — Oosawa size distribution (simplified)

$$\frac{dA_s}{dt} = k_n M^{n_n} \;-\; \underbrace{\frac{k_g M}{n_s - n_n}}_{\text{size-crossing rate}} \cdot A_s$$

$$\frac{dA_i}{dt} = \frac{k_g M}{n_s - n_n} \cdot A_s \;-\; \frac{k_g M}{n_i - n_s} \cdot A_i$$

The term $k_g M / (n_s - n_n)$ is the rate at which each soluble aggregate crosses the soluble/insoluble size threshold by adding monomers (Oosawa elongation). When PQC fails and $M$ rises, this crossing rate accelerates — As converts to Ai faster — creating the **spike then decline** in soluble aggregates.

### Parameter table

| Parameter | ODE value | Meaning |
|-----------|-----------|---------|
| `p` | 400 | Monomer production rate |
| `cT` | 2.0 | Baseline PQC clearance capacity |
| `K` | 17.0 | EC50 for Ai-mediated PQC inhibition (Hill threshold) |
| `h` | 3 | Hill cooperativity of inhibition |
| `kn` | 3 × 10⁻¹³ | Nucleation rate constant |
| `nn` | 6 | Nucleus size (Oosawa / Finke-Watzky) |
| `kg` | 0.05 | Oosawa elongation rate |
| `ns` | 56 | Soluble/insoluble size threshold (`ns − nn = 50` steps) |
| `ni` | 256 | Max insoluble size (`ni − ns = 200` steps) |

Time in **weeks**; concentrations in normalised units matching the Krohn paper (ng/mg protein) in shape.

---

## Gillespie SSA version

The stochastic Gillespie version uses `nn = 2` (tractable nucleus size) with `kn` re-scaled to preserve the same timing as the `nn = 6` ODE. The five reactions are:

| Reaction | Propensity | Biology |
|----------|-----------|---------|
| ∅ → M | `p` | Monomer production |
| M → ∅ | `[cT/(1+(Ai/K)^h)] · M` | **PQC clearance — positive feedback here** |
| 2M → As | `kn · M(M−1)/2` | Nucleation (slow, stochastic) |
| M + As → Ai | `[kg·M/(ns−nn)] · As` | Oosawa growth: As crosses size threshold |
| Ai → ∅ | `[kg·M/(ni−ns)] · Ai` | Ai grows into plaque (leaves tracked pool) |

The stochasticity of nucleation (R3) naturally introduces **run-to-run variability in the spike timing** — a feature absent from the ODE but biologically realistic (different cells/animals show onset at different ages).

---

## Figures produced

| Figure | Description |
|--------|-------------|
| **Fig 1** | Main result: ODE soluble As (left) and insoluble Ai (right) for base, +11% cT, and no-feedback. Reproduces Krohn fig. 4 qualitatively. |
| **Fig 2** | Mechanistic breakdown: PQC efficiency vs time, monomer pool M(t), and phase portrait As vs Ai showing the positive-feedback spiral. |
| **Fig 3** | Gillespie SSA: 12 stochastic runs per scenario (base, +11%, no feedback) for both As and Ai. Shows run-to-run variability in spike timing. |
| **Fig 4** | Sensitivity scan: spike timing and final Ai burden vs cT (85%–130%). Reveals the **non-linear tipping-point sensitivity** — the key argument for targeting cT as intervention. |

---

## Key results

```
ODE model:
  Base  (ABCC1−/−)         spike at t = 18.9 wk    Ai(35 wk) = 43.3
  cT +11%  (ABCB1−/−)      spike at t = 22.8 wk    Ai(35 wk) = 41.4
  No positive feedback      no spike                Ai(35 wk) = 22.6

Effect of +11% cT:  +3.9 wk delay  (21% relative delay in onset)
```

---

## Answering the week question: where to intervene?

Fig. 4 (sensitivity scan) gives the clearest answer. The spike timing as a function of cT is **highly non-linear near the baseline**: a small increase in clearance produces a disproportionately large delay in disease onset. This is the signature of a system near a tipping point, and it is precisely what makes the PQC clearance pathway such an attractive drug target.

The four intervention points from the model, ordered by the leverage they provide:

**1. Increase `cT` — boost PQC clearance capacity.**  
The primary target. In the Krohn model, `cT` aggregates the contributions of α-secretase, autophagy, chaperones, proteasome, LRP1/RAGE, and ABC transporters (ABCC1 being the most effective). The drug thiethylperazine activates ABCC1 transport activity by ~69% in vitro; in the mouse model it produces a significant reduction in soluble and insoluble Aβ. Even an 11% improvement in total `cT` shifts the spike by several weeks because the system is operating near the tipping point.

**2. Reduce `K` — prevent PQC from being disabled by insoluble aggregates.**  
`K` is the EC50 for Ai-mediated inhibition of PQC. Lowering `K` means PQC remains efficient even as insoluble aggregates accumulate, delaying or preventing the runaway positive feedback. Biologically this might correspond to compounds that prevent proteasomal sequestration by fibrils.

**3. Decrease `kn` — inhibit nucleation.**  
The nucleation term `kn · M^nn` is the initiating step. Because `nn = 6`, small reductions in M or kn have a large (sixth-power) effect on the nucleation flux. Anti-nucleation compounds (e.g. EGCG, molecular tweezers) act here.

**4. Decrease `kg` — slow down Oosawa elongation.**  
Slower growth means the size-crossing rate `kg·M/(ns−nn)` is smaller, As stays soluble longer, and the feedback loop through Ai takes longer to build. Fibril capping agents act here.

The critical insight from both the model and the Krohn paper is that **(1) and (2) are downstream of the positive feedback node**, meaning their effect is amplified. Targeting nucleation (3) or elongation (4) acts upstream but the dose–response is more linear.

---

## Connection to lecture material

| Concept | Implementation |
|---------|----------------|
| Positive feedback: Ai --|PQC | Hill function `cT/(1+(Ai/K)^h)` in dM/dt |
| Finke-Watzky 2-step (Week 1) | Nucleation `kn·M^nn` + growth `kg·M·As` |
| Oosawa size distribution | Size-crossing rate `kg·M/(ns−nn)` separates As and Ai pools |
| Krohn 2-module model | Module 1 = dM/dt equation; Module 2 = dAs/dt, dAi/dt |
| 11% cT prediction (Krohn fig 4C) | `run_ode({'cT': cT * 1.11})` reproduces rightward shift |
| Sneppen proteasome oscillations | Same positive feedback logic; here Ai saturates PQC, there fibrils saturate proteasome |
| Separation of timescales | As reaches quasi-SS fast; Ai grows slowly; feedback kicks in on intermediate timescale |

---

## References

- Krohn, M. et al. (2011). *Cerebral amyloid-β proteostasis is regulated by the membrane transport protein ABCC1 in mice.* **Journal of Clinical Investigation**, 121(10), 3924–3931. doi:10.1172/JCI57867
- Sneppen, K. et al. (2009). *Modeling proteasome dynamics in Parkinson's disease.* **Physical Biology**, 6, 036005.
- Lecture slides: *Physics of Molecular Diseases – Week 2, Protein Quality Control* (Prof. Ala Trusina, Niels Bohr Institute, 2020).
- Finke & Watzky (1997). Two-step nucleation-growth model — adapted for protein aggregation throughout this course.
