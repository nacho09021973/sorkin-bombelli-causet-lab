# Gamma resonance hypotheses

Status: exploratory hypothesis note, not a result.

## Context

In SORKIN-2, gamma=0.8 appeared useful in N=12 recovery diagnostics. This note does not claim gamma=0.8 is a physical constant. The aim is to collect possible algorithmic interpretations that may guide future targeted tests.

## Multiplicative cooling

The annealing schedule uses multiplicative cooling:

```text
T_{k+1} = gamma T_k
```

Here gamma controls thermal memory from one step to the next. Larger gamma retains more of the previous temperature, while smaller gamma dissipates thermal mobility more quickly.

## Thermal half-life heuristic

- 0.8^3 = 0.512, close to 1/2.
- 2^(-1/3) = 0.7937005..., close to 0.8.
- Therefore gamma=0.8 approximately corresponds to halving the temperature every three annealing steps.
- This may preserve thermal mobility across a few update layers before freezing.

## Logarithmic cooling rate

- -log(0.8) = 0.223143...
- This is close to a quarter-unit logarithmic decay scale.
- This may be interpreted as a moderate dissipation rate, neither abrupt nor nearly static.

## Accessibility interpretation

- Too small gamma may freeze prematurely.
- Too large gamma may diffuse without committing.
- A middle window near gamma≈0.8 may couple well to the combinatorial rigidity of small N=12 causal orders.
- The relevant object may be the algorithmic accessibility landscape, not causal-set embeddability.

## Future minimal tests

- Compare gamma values around 0.7, 0.8, 0.9 on the same N=12 cases.
- Check whether success is a broad window or a narrow peak.
- Compare topology families: layered, sparse, hub, and sprinkling-like N=12 examples.
- Track whether the same gamma window appears across different seeds.

## Guardrails

- This note is not evidence.
- This note is not a theorem.
- This note does not establish a causal-set invariant.
- Do not cite gamma=0.8 as physically privileged.
- Promote only tested claims to docs/.
- Any future claim must reference named runs or benchmark files.
