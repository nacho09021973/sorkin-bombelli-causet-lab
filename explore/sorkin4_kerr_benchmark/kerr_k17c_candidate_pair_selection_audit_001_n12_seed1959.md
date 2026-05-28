# S4-KERR-K17c candidate-pair selection audit

1. Did K17 select pairs that fail cheap null-admissibility filters?
Yes: many forward-time exterior pairs fail radial/angle gates (named rejections total=168).

2. Does N=12 / seed=1959 contain any plausible near-null candidate pairs?
Selection candidates found=30 across radial-like, low-winding, and sector-aware buckets.

3. Is K18 likely to be informative if run naively?
Probably low-yield without pre-selection. Recommendation: Use radial_time_admissible + low-winding/sector-aware pre-selection before K18; naive K18 is likely low-yield.

4. What selection heuristic should K18 use?
Use deterministic pre-selection: time order + exteriority + radial_time_admissible + low-winding first, then sector-aware fallback.

selection_candidate is not reachability.
rejected_by_selection is not proof of spacelike separation.
residual_probe_pass is not causal reachability.
no production classifier.
no physical/global causal claim.
candidate_undecided remains conservative.
