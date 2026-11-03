# Release tools

`python tools/verify_release.py` checks the release's source/input hashes.
Python is optional for running the model; it is needed only for this integrity
tool or the optional progress GUI. Generated output is intentionally outside
the release manifest and ignored by Git.
