# Pure-`dat` RL shape data

`initial_circle.dat` and `target_naca0012.dat` are the checked-in evaluation
endpoints. The current code reconstructs the circle from the target x-grid,
while retaining both files as explicit teaching data.

Training never modifies these inputs. Checkpoints, histories, rollout `dat`
files, plots, and animations are written to `../output/`.
