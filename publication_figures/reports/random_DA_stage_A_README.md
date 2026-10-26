# Approved Stage A: six random DA tasks, seed 2001

Preflight job: 11137402. Serial pilot array: 11137405, tasks 0–5, global concurrency one. Numerical comparison job: 11137408, dependent on array completion. Report-only job: 11137411, dependent on successful comparison. These dependencies do not submit further DA runs or seeds.

Inputs, manifests, Seurat/BPCells checks and smoke-test outputs: `publication_figures/benchmarks/random_DA_stage_A_2001/preflight/` and `memberships/`.

Random DA outputs: `publication_figures/benchmarks/random_DA_stage_A_2001/runs/{atlas}/seed2001/{02pct,10pct}/`. Each run records the exact command, production-script and membership hashes, job ID and process return code in `execution.json`; R messages are in `DA.log`.

Comparison report and tables: `publication_figures/reports/random_DA_stage_A_2001/`. Accounting comes from Slurm task/step records, including failures. A single random seed does not supply uncertainty across random selections. The fixed common cohort intersects both methods and both budgets within each atlas; it must not be redefined independently per budget.

Once the pilot is complete, regenerate only numerical tables/report on a compute node with:

```bash
sbatch publication_figures/compare_stage_A.sbatch
```

This performs no DA, selection, training or figure generation. Successful DA outputs are protected against overwriting. Do not resubmit the pilot to extend seeds: further experiments require a new explicit approval.
