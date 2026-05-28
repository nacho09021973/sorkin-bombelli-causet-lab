# S4-KERR-K17d cloud-size / seed / spin scan

1. Does best_residual decrease with N?
monotone_N_decrease_flag = No. Median by N: N=12: 4.5421e-01, N=24: 2.6508e-01, N=48: 2.6508e-01.

2. Are there near-hits (W_TOL < residual <= 10*W_TOL)?
near_hits_overall = 0. residual_le_w_tol_overall = 0.

3. Does any cell cross W_TOL or 10*W_TOL?
N_at_which_W_TOL_first_crossed = None.
N_at_which_10x_W_TOL_first_crossed = None.

4. Is K18 informative naively on this generator?
Not informative naively; random cloud-cloud strategy fails to shrink residual with N.

5. Recommendation:
Random cloud-cloud pair selection is not informative on this generator at this scale; consider inspecting the cloud distribution (Option B) or designing structured candidate pairs (Option C) before K18.

6. Per-cell breakdown (rejection counts + best residual + near-hits):

| N | seed | spin | n_sel | rej_time | rej_radial | rej_angular | unresolved | best_residual | near_hits |
|---|---|---|---|---|---|---|---|---|---|
| 12 | 1959 | 0.00 | 10 | 66 | 49 | 7 | 0 | 4.9806e-01 | 0 |
| 12 | 1959 | 0.25 | 10 | 66 | 49 | 7 | 0 | 5.0050e-01 | 0 |
| 12 | 1959 | 0.50 | 10 | 66 | 49 | 7 | 0 | 5.0300e-01 | 0 |
| 12 | 1960 | 0.00 | 9 | 66 | 42 | 15 | 0 | 4.3508e-01 | 0 |
| 12 | 1960 | 0.25 | 9 | 66 | 42 | 15 | 0 | 4.3331e-01 | 0 |
| 12 | 1960 | 0.50 | 9 | 66 | 42 | 15 | 0 | 4.3298e-01 | 0 |
| 12 | 1961 | 0.00 | 15 | 66 | 40 | 11 | 0 | 2.4980e-01 | 0 |
| 12 | 1961 | 0.25 | 13 | 66 | 42 | 11 | 0 | 4.5849e-01 | 0 |
| 12 | 1961 | 0.50 | 13 | 66 | 42 | 11 | 0 | 4.5421e-01 | 0 |
| 24 | 1959 | 0.00 | 39 | 276 | 197 | 40 | 0 | 3.4689e-01 | 0 |
| 24 | 1959 | 0.25 | 39 | 276 | 197 | 40 | 0 | 3.4262e-01 | 0 |
| 24 | 1959 | 0.50 | 38 | 276 | 198 | 40 | 0 | 3.3669e-01 | 0 |
| 24 | 1960 | 0.00 | 44 | 276 | 184 | 48 | 0 | 2.6079e-01 | 0 |
| 24 | 1960 | 0.25 | 44 | 276 | 184 | 48 | 0 | 2.6508e-01 | 0 |
| 24 | 1960 | 0.50 | 43 | 276 | 185 | 48 | 0 | 2.6720e-01 | 0 |
| 24 | 1961 | 0.00 | 46 | 276 | 180 | 50 | 0 | 1.3264e-01 | 0 |
| 24 | 1961 | 0.25 | 44 | 276 | 182 | 50 | 0 | 1.2963e-01 | 0 |
| 24 | 1961 | 0.50 | 44 | 276 | 184 | 48 | 0 | 1.2547e-01 | 0 |
| 48 | 1959 | 0.00 | 154 | 1128 | 824 | 150 | 0 | 2.8346e-01 | 0 |
| 48 | 1959 | 0.25 | 152 | 1128 | 826 | 150 | 0 | 2.8147e-01 | 0 |
| 48 | 1959 | 0.50 | 150 | 1128 | 833 | 145 | 0 | 2.7959e-01 | 0 |
| 48 | 1960 | 0.00 | 195 | 1128 | 730 | 203 | 0 | 2.6079e-01 | 0 |
| 48 | 1960 | 0.25 | 194 | 1128 | 732 | 202 | 0 | 2.6508e-01 | 0 |
| 48 | 1960 | 0.50 | 191 | 1128 | 737 | 200 | 0 | 2.6720e-01 | 0 |
| 48 | 1961 | 0.00 | 216 | 1128 | 717 | 195 | 0 | 1.3264e-01 | 0 |
| 48 | 1961 | 0.25 | 214 | 1128 | 719 | 195 | 0 | 1.2963e-01 | 0 |
| 48 | 1961 | 0.50 | 214 | 1128 | 722 | 192 | 0 | 1.2547e-01 | 0 |

Refinement boundary: K17d does not refine top cases; refinement is reserved for a separate phase.

near_hit is not reachability.
selection_candidate is not reachability.
rejected_by_selection is not proof of spacelike separation.
residual_probe_pass is not causal reachability.
candidate_miss is not proof of spacelike separation.
no production classifier introduced.
no physical/global causal claim introduced.
no Level-B Hawking/Bekenstein claim introduced.
no causal_true/false relations decided.
no production sprinkling touched.
