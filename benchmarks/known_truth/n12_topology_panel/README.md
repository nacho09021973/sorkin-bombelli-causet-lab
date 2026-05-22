# SORKIN-2 N=12 topology panel inputs

This directory contains explicitly constructed N=12 target causal orders for
the SORKIN-2 topology panel.

These files are input artifacts only. They are not annealer results, do not
imply that any run has been executed, and must not be cited as recovery output.

Purpose:

- study how target-poset topology affects historical annealer accessibility;
- provide controlled N=12 anchors before selecting additional random
  sprinklings;
- keep topology and recoverability questions separate.

Cases:

- `chain_12_d2`: total chain on 12 elements; dense relation-count extreme.
- `antichain_12_d2`: 12 unrelated elements; sparse all-spacelike extreme.
- `layered_4_4_4_d2`: three antichain layers of four elements each, with full
  order between layers; controlled hub/layer-rich topology.

These cases do not prove manifoldlikeness, do not test general embeddability,
and do not prove completeness of the Bombelli annealer. Future runs should
treat them as constructed-truth or structural accessibility diagnostics.
