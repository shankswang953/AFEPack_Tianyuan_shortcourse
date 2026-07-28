#include "airfoil_bezier.hpp"

#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;
using airfoil_demo::AirfoilData;
using airfoil_demo::FittedAirfoil;

int main(int argc, char** argv)
{
  try {
    if (argc < 3 || argc > 9) {
      std::cerr
          << "usage: generate_airfoil_geometry DATA_FILE OUTPUT_DIR"
          << " [boundary_points=96] [motion_center=0.35]"
          << " [motion_width=0.14] [upper_shift=0.030]"
          << " [lower_shift=0.010] [MOVED_DATA_FILE]\n"
          << "If MOVED_DATA_FILE is supplied, it defines the target"
          << " airfoil and the Gaussian parameters are ignored.\n";
      return 1;
    }

    const fs::path data_file = argv[1];
    const fs::path output_dir = argv[2];
    const int boundary_points = argc > 3 ? std::stoi(argv[3]) : 96;
    const double motion_center = argc > 4 ? std::stod(argv[4]) : 0.35;
    const double motion_width = argc > 5 ? std::stod(argv[5]) : 0.14;
    const double upper_shift = argc > 6 ? std::stod(argv[6]) : 0.030;
    const double lower_shift = argc > 7 ? std::stod(argv[7]) : 0.010;
    const bool use_moved_data_file = argc > 8;
    const fs::path moved_data_file =
        use_moved_data_file ? fs::path(argv[8]) : fs::path();

    constexpr int control_points_per_surface = 15;
    // The training case starts from a circle.  A weak regularization retains
    // that circle while still damping tiny accumulated data perturbations.
    constexpr double fitting_smoothness = 0.01;
    constexpr double outer_center_x = 0.5;
    constexpr double outer_center_y = 0.0;
    constexpr double outer_radius = 4.0;
    constexpr int outer_points = 64;

    fs::create_directories(output_dir);

    AirfoilData initial =
        airfoil_demo::read_uiuc_airfoil(data_file);
    airfoil_demo::close_airfoil_endpoints(initial);
    AirfoilData moved;
    if (use_moved_data_file) {
      moved = airfoil_demo::read_uiuc_airfoil(moved_data_file);
      airfoil_demo::close_airfoil_endpoints(moved);
      if (moved.upper.size() != initial.upper.size()
          || moved.lower.size() != initial.lower.size()) {
        throw std::runtime_error(
            "initial and moved airfoil data have different sizes");
      }
      for (std::size_t i = 0; i < initial.upper.size(); ++i) {
        if (std::abs(initial.upper[i].x - moved.upper[i].x) > 1.0e-12) {
          throw std::runtime_error(
              "upper-surface x coordinates changed in MOVED_DATA_FILE");
        }
      }
      for (std::size_t i = 0; i < initial.lower.size(); ++i) {
        if (std::abs(initial.lower[i].x - moved.lower[i].x) > 1.0e-12) {
          throw std::runtime_error(
              "lower-surface x coordinates changed in MOVED_DATA_FILE");
        }
      }
    } else {
      moved = airfoil_demo::apply_gaussian_motion(
          initial,
          motion_center,
          motion_width,
          upper_shift,
          lower_shift);
    }

    const FittedAirfoil initial_fit = airfoil_demo::fit_and_sample(
        initial,
        control_points_per_surface,
        fitting_smoothness,
        boundary_points);
    const FittedAirfoil moved_fit = airfoil_demo::fit_and_sample(
        moved,
        control_points_per_surface,
        fitting_smoothness,
        boundary_points);

    airfoil_demo::write_easymesh_domain(
        output_dir / "airfoil.d",
        initial_fit.boundary_points,
        outer_center_x,
        outer_center_y,
        outer_radius,
        outer_points);
    airfoil_demo::write_indexed_points(
        output_dir / "boundary_initial.dat",
        initial_fit.boundary_points);
    airfoil_demo::write_indexed_points(
        output_dir / "boundary_moved.dat",
        moved_fit.boundary_points);
    airfoil_demo::write_fit_csv(
        output_dir / "fit_initial.csv",
        initial,
        initial_fit);
    airfoil_demo::write_fit_csv(
        output_dir / "fit_moved.csv",
        moved,
        moved_fit);
    airfoil_demo::write_uiuc_airfoil(
        output_dir / "airfoil_moved.dat",
        moved);

    double maximum_boundary_motion = 0.0;
    for (std::size_t i = 0; i < initial_fit.boundary_points.size(); ++i) {
      maximum_boundary_motion = std::max(
          maximum_boundary_motion,
          airfoil_demo::distance(
              initial_fit.boundary_points[i],
              moved_fit.boundary_points[i]));
    }

    std::cout << "Input surface points: "
              << initial.upper.size() << " upper + "
              << initial.lower.size() << " lower\n";
    std::cout << "Bezier control points: "
              << control_points_per_surface << " per surface\n";
    std::cout << "EasyMesh airfoil points: "
              << initial_fit.boundary_points.size() << "\n";
    std::cout << "Maximum prescribed boundary motion: "
              << maximum_boundary_motion << "\n";
    std::cout << "Target geometry: "
              << (use_moved_data_file ? moved_data_file.string()
                                      : "built-in Gaussian test motion")
              << "\n";
    std::cout << "EasyMesh input: " << (output_dir / "airfoil.d") << "\n";
  } catch (const std::exception& error) {
    std::cerr << "generate_airfoil_geometry: " << error.what() << "\n";
    return 2;
  }
  return 0;
}
