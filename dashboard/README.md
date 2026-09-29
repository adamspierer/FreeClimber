# dashboard/

Lightweight exploratory scripts for FreeClimber-FNG outputs. These are
helpers, not a finished UI — they wrap the existing CSV outputs in quick
matplotlib summaries.

## Scripts

- `plot_tortuosity.py` — Histograms of the three per-bout tortuosity
  metrics from one or more `*.tortuosity.csv` files. Headless (Agg backend);
  safe to run in a non-interactive shell.

  ```
  python dashboard/plot_tortuosity.py path/to/<video>.tortuosity.csv
  ```

  Writes `tortuosity_summary.png` next to the first input CSV (or to
  `--out` if specified).

The per-fly tortuosity outputs are produced only when
`analysis_mode='individual'` is set in the configuration file; see the
*Per-fly tortuosity metrics* section of the top-level `README.md`.
