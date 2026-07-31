# README result figures

This directory contains selected high-resolution static results embedded in
the numbered example README files. These are presentation copies of
reproducible files below the examples' ignored `output/` directories.

| example | tracked figure | generated source |
|---|---|---|
| `00` | `00_poisson_solution.png` | `00_poisson_basic/output/u.dx`, rendered with `common/mesh/visualize_dx.py` |
| `01` | `01_indicator_scaling.png` | `01_local_refinement/output/figures/01_indicator_scaling.png` |
| `01` | `01_refinement_comparison.png` | `01_local_refinement/output/figures/02_refinement_comparison.png` |
| `02` | `02_tianyuan_merge_overview.png` | `02_mesh_merge/output/figures/tianyuan_merge_overview.png` |
| `05` | `05_final_meshes_triptych.png` | the 10/10/7 comparison's `teaching/final_meshes_triptych.png` |
| `05` | `05_functional_error_vs_dofs.png` | the 10/10/7 comparison's `functional_error_vs_dofs.png` |
| `06` | `06_mesh_motion_overview.png` | `06_airfoil_mesh_motion/output/figures/mesh_motion_overview.png` |
| `07` | `07_shape_history.svg` | `07_poisson_shape_optimization/output/shape_history.svg` |
| `08` | `08_ml_iteration_final_mesh.png` | `08_airfoil_ml_iteration/output/policy_rollout/boundary_mesh_evolution_final.png` |
| `09` | `09_fixed_topology_mesh_path.png` | `09_airfoil_barycentric_motion/output/figures/02_fixed_topology_mesh_path.png` |
| `09` | `09_mesh_quality_along_path.png` | `09_airfoil_barycentric_motion/output/figures/03_mesh_quality_along_path.png` |
| `10` | `10_training_history.png` | `10_airfoil_rl_dat/output/training_history.png` |
| `10` | `10_evaluation_loss.png` | `10_airfoil_rl_dat/output/evaluation/loss_history.png` |

Examples `03` and `04` already keep their primary figures in tracked local
`figures/` directories, so their README files reference those originals
instead of duplicating them here.
