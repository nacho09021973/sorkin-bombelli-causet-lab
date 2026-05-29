# S4-KERR-K18B single-candidate geometric audit 001

## Status

Completed (local diagnostic audit only, one K18A family only).

## Input artifacts actually used

- `explore/sorkin4_kerr_benchmark/kerr_k18a_topk_candidate_sandbox_001.md`
- `explore/sorkin4_kerr_benchmark/kerr_k17d_cloud_size_seed_scan_001.csv`
- `explore/sorkin4_kerr_benchmark/kerr_k17e_necessary_causal_filter_audit_001.csv`
- `explore/sorkin4_kerr_benchmark/audit_kerr_k17d_cloud_size_seed_scan_001.py`
- `explore/sorkin4_kerr_benchmark/kerr_k17d_cloud_size_seed_scan_001.md`

## Candidate family audited

- `N=24`, `seed=1961`, `event_A=15`, `event_B=4`
- `spin_a in {0.00, 0.25, 0.50}`
- K18A best-family rows:
  - `k17d_N24_seed1961_a0.00_A15_B4`
  - `k17d_N24_seed1961_a0.25_A15_B4`
  - `k17d_N24_seed1961_a0.50_A15_B4`

## Exact commands used

```bash
git status --short
```

```bash
rg --files explore/sorkin4_kerr_benchmark
```

```bash
rg -n "N=24|seed=1961|event_A=15|event_B=4|0.125469416355373|k18a" \
  explore/sorkin4_kerr_benchmark/kerr_k18a_topk_candidate_sandbox_001.md \
  explore/sorkin4_kerr_benchmark/kerr_k17e_necessary_causal_filter_audit_001.csv \
  explore/sorkin4_kerr_benchmark/kerr_k17e_necessary_causal_filter_audit_001.md \
  explore/sorkin4_kerr_benchmark/kerr_k17d_cloud_size_seed_scan_001.csv
```

```bash
sed -n '1,260p' explore/sorkin4_kerr_benchmark/kerr_k18a_topk_candidate_sandbox_001.md
```

```bash
sed -n '1,220p' explore/sorkin4_kerr_benchmark/audit_kerr_k17d_cloud_size_seed_scan_001.py
```

```bash
sed -n '1,220p' explore/sorkin4_kerr_benchmark/kerr_k17d_cloud_size_seed_scan_001.md
```

```bash
python3 - <<'PY'
import csv
from pathlib import Path
p=Path('explore/sorkin4_kerr_benchmark/kerr_k17d_cloud_size_seed_scan_001.csv')
rows=[]
with p.open() as f:
    for r in csv.DictReader(f):
        if int(r['N'])==24 and int(r['seed'])==1961 and int(r['event_A_index'])==15 and int(r['event_B_index'])==4 and float(r['spin_a']) in (0.0,0.25,0.5):
            rows.append(r)
for x in sorted(rows,key=lambda z:float(z['spin_a'])):
    print(x['case_id'],x['best_residual'],x['best_b'],x['best_lambda'],x['best_direction'],x['best_sector_m'])
print('count',len(rows))
base=sorted(rows,key=lambda z:float(z['spin_a']))
if len(base)==3:
    r0=float(base[0]['best_residual']); r1=float(base[1]['best_residual']); r2=float(base[2]['best_residual'])
    print('delta(0.00->0.25)=',r1-r0)
    print('delta(0.25->0.50)=',r2-r1)
    print('delta(0.00->0.50)=',r2-r0)
PY
```

```bash
python3 - <<'PY'
import csv
from pathlib import Path
p=Path('explore/sorkin4_kerr_benchmark/kerr_k17e_necessary_causal_filter_audit_001.csv')
rows=[]
with p.open() as f:
    for r in csv.DictReader(f):
        if int(r['N'])==24 and int(r['seed'])==1961 and int(r['event_A_index'])==15 and int(r['event_B_index'])==4 and float(r['spin_a']) in (0.0,0.25,0.5):
            rows.append(r)
for x in sorted(rows,key=lambda z:float(z['spin_a'])):
    print(x['case_id'],x['best_residual_k17d'],x['combined_not_excluded'])
print('count',len(rows))
PY
```

## Observed residual pattern across spin_a

For this fixed candidate family (`N=24`, `seed=1961`, `A15->B4`), residual decreases monotonically as `spin_a` increases:

- `a=0.00`: `best_residual=0.13263756990470288`
- `a=0.25`: `best_residual=0.12963041412301557` (delta vs 0.00: `-0.0030071557816873096`)
- `a=0.50`: `best_residual=0.12546941635537312` (delta vs 0.25: `-0.004160997767642449`)

Net drift over `a=0.00 -> 0.50`: `-0.007168153549329759`.

Observed shape is a smooth monotone drift across the sampled spin grid, not a discontinuous jump.

## Available geometric/provenance data

From K17d CSV for the audited rows:

- Shared provenance:
  - `N=24`, `seed=1961`, `event_A_index=15`, `event_B_index=4`
  - `delta_t_AB=0.9879731090347068`
  - `angular_separation_mod_2pi=0.14930012223885214`
  - `best_direction=ingoing`, `best_sector_m=0`
- Coordinate trend with spin:
  - `t_A` and `t_B` unchanged across spin (`2.9103567417252916`, `3.8983298507599984`)
  - `phi_A` and `phi_B` unchanged across spin (`4.9572224569140015`, `5.106522579152854`)
  - `r_A` decreases with spin: `4.931311593949162 -> 4.927942433808879 -> 4.917680370166319`
  - `r_B` decreases with spin: `4.570525525651426 -> 4.565603803055676 -> 4.550588789137428`
  - `delta_r_AB` becomes slightly more negative: `-0.3607860682977355 -> -0.3623386307532037 -> -0.36709158102889106`

From K17e CSV for the same rows:

- `combined_not_excluded=True` for all three spin values.
- K17e preserves the same K17d best residual values for this family.

Grid-position check from existing K17d script (`audit_kerr_k17d_cloud_size_seed_scan_001.py`):

- `PROBE_B_GRID = (-1.0, -0.5, 0.0, 0.5, 1.0)`
- `PROBE_LAMBDA_GRID = (0.5, 1.0, 2.0)`
- audited best values are `best_b=1.0` and `best_lambda=0.5` for all three spins.

Therefore, in the K17d probe grid:

- `best_b=1.0` is a boundary value (upper edge).
- `best_lambda=0.5` is a boundary value (lower edge).

## Interpretation

Within this single family, low residual is consistent across `spin_a={0.00,0.25,0.50}` with smooth monotone drift, and the selected probe controls (`best_b`, `best_lambda`) remain fixed at grid boundaries. This supports a local statement: this family behaves as a stable low-residual pocket under the sampled spin sweep, with boundary-attained probe parameters.

No further inference is made beyond this local geometric/probe behavior.

## Guardrails

- audited candidate != causal_true
- low residual != reachability
- residual stability across spin_a != physical proof
- candidate miss/hit != spacelike/timelike classification
- this audit is local to one K18A family

## Verdict

K18B local audit completed for the requested family only. The low residual appears spin-stable in a smooth drift sense across `{0.00,0.25,0.50}`, while `best_b=1.0` and `best_lambda=0.5` are boundary grid selections in the existing K17d probe lattice.

Because both probe controls are boundary-attained, this does not justify a positive K18C causal/geometric claim. At most, it motivates a local boundary-sensitivity check if this branch is continued.
