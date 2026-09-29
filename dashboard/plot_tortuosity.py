#!/usr/bin/env python
# -*- coding: utf-8 -*-

## File name : plot_tortuosity.py
## Purpose   : Lightweight exploratory plots for *.tortuosity.csv outputs of
##             FreeClimber-FNG individual-mode runs.

"""
Plot per-bout tortuosity metrics for one or more *.tortuosity.csv files.

Usage
    python dashboard/plot_tortuosity.py <tortuosity_csv> [<tortuosity_csv> ...]
        [--out <output_png>]

If no --out is supplied, a PNG named tortuosity_summary.png is written next
to the first input CSV.

This script is intentionally minimal -- it is a starting point for the
analysis dashboard, not a finished UI. It does not modify FreeClimber's own
outputs.

The CSV schema expected is the one written by detector.compute_tortuosity():
    vial, event_idx, particle, frame_start, frame_end, n_points,
    path_length, net_displacement, tortuosity, straightness,
    mean_turning_angle_rad
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # non-interactive backend; required for headless runs
import matplotlib.pyplot as plt


def _load(paths):
    """Concatenate one-or-more tortuosity CSVs into a single labeled frame."""
    frames = []
    for p in paths:
        df = pd.read_csv(p)
        df['_source'] = os.path.basename(p)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def _safe_auto_bins(values, n=30):
    """Build a bin spec robust to degenerate inputs (empty / single-value).

    Matplotlib's automatic bin sizing fails when min == max; widen the range
    by a small symmetric margin in that case.
    """
    if len(values) == 0:
        return np.linspace(0.0, 1.0, n + 1)
    lo, hi = float(np.min(values)), float(np.max(values))
    if lo == hi:
        pad = 0.5 if lo == 0 else abs(lo) * 0.05
        return np.linspace(lo - pad, hi + pad, n + 1)
    return n


def _plot(df, out_path):
    """Render a 3-panel summary: tortuosity, straightness, mean turning angle."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))

    # Drop rows where the metric is undefined (closed bouts produce NaN T).
    tort = df['tortuosity'].dropna()
    straight = df['straightness'].dropna()
    theta = df['mean_turning_angle_rad'].dropna()

    axes[0].hist(tort, bins=_safe_auto_bins(tort),
                 color='#4477aa', edgecolor='white')
    axes[0].set_xlabel('classical tortuosity (L / D)')
    axes[0].set_ylabel('bout-particle count')
    axes[0].set_title('Tortuosity (n=%d)' % len(tort))

    axes[1].hist(straight, bins=np.linspace(0.0, 1.0, 31),
                 color='#228833', edgecolor='white')
    axes[1].set_xlabel('straightness (D / L)')
    axes[1].set_title('Straightness (n=%d)' % len(straight))

    axes[2].hist(theta, bins=np.linspace(0.0, np.pi, 31),
                 color='#ee6677', edgecolor='white')
    axes[2].set_xlabel('mean |turning angle| (rad)')
    axes[2].set_title('Turning angle (n=%d)' % len(theta))

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print('Saved: %s' % out_path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument('csv', nargs='+', help='one or more *.tortuosity.csv files')
    parser.add_argument('--out', default=None,
                        help='output PNG path (default: next to the first input)')
    args = parser.parse_args(argv)

    df = _load(args.csv)
    if df.empty:
        print('No rows in any input CSV; nothing to plot.', file=sys.stderr)
        return 1

    out = args.out or os.path.join(os.path.dirname(os.path.abspath(args.csv[0])),
                                   'tortuosity_summary.png')
    _plot(df, out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
