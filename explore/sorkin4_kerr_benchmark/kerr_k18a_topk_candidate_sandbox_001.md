# S4-KERR-K18A top-K candidate sandbox

Status: diagnostic triage artifact (not a causal classifier, not a new filter).

## Inputs found and used

- `explore/sorkin4_kerr_benchmark/kerr_k17d_cloud_size_seed_scan_001.csv`
- `explore/sorkin4_kerr_benchmark/kerr_k17e_necessary_causal_filter_audit_001.csv`
- Context only (not parsed for table values):
  - `explore/sorkin4_kerr_benchmark/kerr_k17d_cloud_size_seed_scan_001.md`
  - `explore/sorkin4_kerr_benchmark/kerr_k17e_necessary_causal_filter_audit_001.md`

## Exact command used

```bash
python3 - <<'PY'
import csv
from pathlib import Path

k17d=Path('explore/sorkin4_kerr_benchmark/kerr_k17d_cloud_size_seed_scan_001.csv')
k17e=Path('explore/sorkin4_kerr_benchmark/kerr_k17e_necessary_causal_filter_audit_001.csv')

rows=[]
with k17d.open() as f:
    for r in csv.DictReader(f):
        r['best_residual']=float(r['best_residual'])
        r['N']=int(r['N'])
        r['seed']=int(r['seed'])
        r['spin_a']=float(r['spin_a'])
        rows.append(r)

survivor_case_ids=set()
with k17e.open() as f:
    for r in csv.DictReader(f):
        if r['combined_not_excluded'].strip().lower() in ('true','1','yes'):
            survivor_case_ids.add(r['case_id'])

rows=[r for r in rows if r['case_id'] in survivor_case_ids]

def dedup_signature(r):
    return (
        r['seed'], f"{r['spin_a']:.2f}",
        f"{float(r['t_A']):.9f}", f"{float(r['r_A']):.9f}", f"{float(r['phi_A']):.9f}",
        f"{float(r['t_B']):.9f}", f"{float(r['r_B']):.9f}", f"{float(r['phi_B']):.9f}",
        f"{r['best_residual']:.12f}",
        f"{float(r['best_b']):.6f}", f"{float(r['best_lambda']):.6f}",
        r['best_direction'], str(int(float(r['best_sector_m'])))
    )

best_for_signature={}
for r in rows:
    s=dedup_signature(r)
    cur=best_for_signature.get(s)
    if cur is None or (r['best_residual'], r['N'], r['case_id']) < (cur['best_residual'], cur['N'], cur['case_id']):
        best_for_signature[s]=r

dedup_sorted=sorted(
    best_for_signature.values(),
    key=lambda r:(r['best_residual'], r['N'], r['seed'], r['spin_a'], r['case_id'])
)

print('rows_in=',len(rows))
print('rows_after_dedup=',len(dedup_sorted))
for i,r in enumerate(dedup_sorted[:25],1):
    print(i, r['case_id'], r['N'], r['seed'], r['spin_a'], f"{r['best_residual']:.15f}")
PY
```

## Deduplication rule

This K18A artifact treats `N=24` and `N=48` as candidate sources, not independent confirmation.
A row is considered an obvious clone if it shares the same geometric/probe signature after deterministic numeric rounding and differs only by sampling/source details (including `N` and event index labels).

Signature fields:

- `seed`, `spin_a`
- endpoint coordinates: `t_A`, `r_A`, `phi_A`, `t_B`, `r_B`, `phi_B` (rounded to 9 decimals)
- K17d probe result fields: `best_residual` (12 decimals), `best_b`, `best_lambda`, `best_direction`, `best_sector_m`

Tie-break inside each clone group: keep the row with smallest tuple `(best_residual, N, case_id)`.

Observed counts from this rule:

- K17d rows matched to K17E `combined_not_excluded=true`: 2159
- deduplicated rows: 1680
- obvious clones removed: 479

## Deduplicated top-K candidates (K=25)

| rank | case_id | N | seed | spin_a | best_residual | event_A | event_B | best_b | best_lambda | best_direction | best_sector_m | source_files |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---|
| 1 | k17d_N24_seed1961_a0.50_A15_B4 | 24 | 1961 | 0.50 | 0.125469416355373 | 15 | 4 | 1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 2 | k17d_N24_seed1961_a0.25_A15_B4 | 24 | 1961 | 0.25 | 0.129630414123016 | 15 | 4 | 1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 3 | k17d_N24_seed1961_a0.00_A15_B4 | 24 | 1961 | 0.00 | 0.132637569904703 | 15 | 4 | 1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 4 | k17d_N24_seed1961_a0.50_A18_B13 | 24 | 1961 | 0.50 | 0.140210249466766 | 18 | 13 | 1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 5 | k17d_N24_seed1961_a0.25_A18_B13 | 24 | 1961 | 0.25 | 0.143704583863432 | 18 | 13 | 1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 6 | k17d_N24_seed1961_a0.00_A18_B13 | 24 | 1961 | 0.00 | 0.146290557740909 | 18 | 13 | -1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 7 | k17d_N24_seed1961_a0.50_A23_B12 | 24 | 1961 | 0.50 | 0.209695069865228 | 23 | 12 | 1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 8 | k17d_N24_seed1961_a0.25_A23_B12 | 24 | 1961 | 0.25 | 0.213511624500655 | 23 | 12 | 1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 9 | k17d_N24_seed1961_a0.00_A23_B12 | 24 | 1961 | 0.00 | 0.216861678053426 | 23 | 12 | -1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 10 | k17d_N48_seed1961_a0.00_A34_B45 | 48 | 1961 | 0.00 | 0.229438958670313 | 34 | 45 | -1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 11 | k17d_N48_seed1961_a0.25_A34_B45 | 48 | 1961 | 0.25 | 0.231920762106059 | 34 | 45 | -1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 12 | k17d_N48_seed1961_a0.50_A34_B45 | 48 | 1961 | 0.50 | 0.234464022559513 | 34 | 45 | -1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 13 | k17d_N12_seed1961_a0.00_A5_B8 | 12 | 1961 | 0.00 | 0.249804133143940 | 5 | 8 | 0.0 | 0.5 | outgoing | 0 | k17d CSV + k17e CSV |
| 14 | k17d_N24_seed1960_a0.00_A9_B15 | 24 | 1960 | 0.00 | 0.260789669779916 | 9 | 15 | -1.0 | 1.0 | ingoing | 0 | k17d CSV + k17e CSV |
| 15 | k17d_N24_seed1960_a0.25_A9_B15 | 24 | 1960 | 0.25 | 0.265076632832531 | 9 | 15 | -1.0 | 1.0 | ingoing | 0 | k17d CSV + k17e CSV |
| 16 | k17d_N24_seed1960_a0.50_A9_B15 | 24 | 1960 | 0.50 | 0.267196125059212 | 9 | 15 | -1.0 | 1.0 | ingoing | 0 | k17d CSV + k17e CSV |
| 17 | k17d_N48_seed1959_a0.50_A12_B25 | 48 | 1959 | 0.50 | 0.279593767770982 | 12 | 25 | 1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 18 | k17d_N48_seed1959_a0.25_A12_B25 | 48 | 1959 | 0.25 | 0.281466969532866 | 12 | 25 | 1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 19 | k17d_N48_seed1959_a0.00_A12_B25 | 48 | 1959 | 0.00 | 0.283458219378816 | 12 | 25 | 1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 20 | k17d_N48_seed1961_a0.50_A9_B46 | 48 | 1961 | 0.50 | 0.286698708013317 | 9 | 46 | 1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 21 | k17d_N48_seed1961_a0.25_A9_B46 | 48 | 1961 | 0.25 | 0.288585457094808 | 9 | 46 | 1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 22 | k17d_N48_seed1961_a0.00_A9_B46 | 48 | 1961 | 0.00 | 0.290788575684682 | 9 | 46 | -1.0 | 0.5 | ingoing | 0 | k17d CSV + k17e CSV |
| 23 | k17d_N24_seed1960_a0.50_A14_B0 | 24 | 1960 | 0.50 | 0.292545302087444 | 14 | 0 | 1.0 | 0.5 | outgoing | 0 | k17d CSV + k17e CSV |
| 24 | k17d_N24_seed1960_a0.25_A14_B0 | 24 | 1960 | 0.25 | 0.295139387011575 | 14 | 0 | 1.0 | 0.5 | outgoing | 0 | k17d CSV + k17e CSV |
| 25 | k17d_N48_seed1961_a0.50_A41_B43 | 48 | 1961 | 0.50 | 0.295315481478803 | 41 | 43 | 1.0 | 1.0 | outgoing | 0 | k17d CSV + k17e CSV |

## Caveats

- top-K candidate != causal_true
- candidate_hit != reachability
- candidate_miss != spacelike
- low residual != proof
- deduplicated clone != independent evidence
- not_excluded is not reachability
- rejected_by_filter is not proof of spacelike separation
- no causal_true / causal_false labels are assigned here
- no production classifier is introduced
- no physical/global causal claim is introduced

## Operational verdict

Some candidates merit K18B geometric audit.

Reasoning constrained to current data: K17E is non-selective (no additional exclusion), while K18A produces a deterministic, provenance-preserving, clone-aware top-K ordering from K17d outputs that can be audited geometrically without treating cloud-cloud i.i.d. sampling as confirmation.
