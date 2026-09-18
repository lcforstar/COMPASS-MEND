# GitHub upload package

This directory is a curated copy of the COMPASS-MEND working tree prepared on
2026-09-18. It contains the current MATLAB production source and the input data
needed by the retained site packages. It is not a Git clone and contains no Git
history.

## Included

- `src/`: current production MATLAB source, including uncommitted working-tree changes.
- Root launchers and `scripts/`: normal run and progress-monitor entry points.
- `data/sites/`: stable Changwu, Changwu-Apple, and laboratory water site packages.
- Shared runtime inputs under `data/atmosphere`, `data/forcing`, `data/global`,
  `data/input`, `data/model_library_nc`, and `data/observations`.
- Lightweight Markdown/SVG/JSON technical documentation.

## Deliberately excluded

- Git history and local agent configuration.
- Generated `output/` and `outputs/` trees, checkpoints, logs, figures, and run archives.
- Tests, verification fixtures, diagnostic/calibration tools, plotting scripts, and
  temporary experimental runners.
- Timestamped calibration/run-specific site copies.
- PDF literature caches, temporary downloads, backup files, and archived plans.

The largest included file is below GitHub's 100 MB per-file limit. Before making
the repository public, confirm that you are authorized to redistribute every
included third-party input dataset and add the intended software/data licenses.
