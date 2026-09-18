# Outer-step transaction layer

`step` owns the model-wide outer-step contract. It validates the step request,
clips the accepted clock interval, builds and commits destination-keyed pending
state, and constructs the step report.

The package contains no canopy, soil, plant, or radiation physics. Numerical
trial/retry control remains in `solver`; physical kernels remain in their
discipline packages; diagnostics and file output remain in `diagnostics` and
`io`.

State updates are committed only after the soil outer step is accepted. A
rejected step must not advance prognostic state, management events, the runtime
clock, diagnostics, checkpoints, or partial output.

