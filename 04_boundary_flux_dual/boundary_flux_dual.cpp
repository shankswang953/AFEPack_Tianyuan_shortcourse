#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <AFEPack/AMGSolver.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/HGeometry.h>
#include <AFEPack/Operator.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/Vector.h>

namespace {

constexpr int kDimension = 2;
constexpr int kBottomBoundary = 1;
constexpr double kConductivity = 1.0;
constexpr double kManufacturedAmplitude = 0.04;

// Heat flux is observed through this part of the cooled bottom boundary.
constexpr double kFluxPatchLeft = 0.40;
constexpr double kFluxPatchRight = 0.60;
constexpr int kUniformRefinementLevel = 2;

double zero_boundary_value(const double*) {
  return 0.0;
}

double exact_temperature(const double* point) {
  const double pi = std::acos(-1.0);
  const double x_factor =
      1.0 - kManufacturedAmplitude * std::cos(2.0 * pi * point[0]);
  return point[1] * (2.0 - point[1]) * x_factor;
}

double manufactured_heat_source(const double* point) {
  const double pi = std::acos(-1.0);
  const double y_factor = point[1] * (2.0 - point[1]);
  return 2.0 - kManufacturedAmplitude *
      (2.0 + 4.0 * pi * pi * y_factor) *
      std::cos(2.0 * pi * point[0]);
}

double exact_flux_functional() {
  const double pi = std::acos(-1.0);
  return 2.0 * (kFluxPatchRight - kFluxPatchLeft) -
      kManufacturedAmplitude / pi *
      (std::sin(2.0 * pi * kFluxPatchRight) -
       std::sin(2.0 * pi * kFluxPatchLeft));
}

double distance_squared(const Point<kDimension>& first,
                        const Point<kDimension>& second) {
  double result = 0.0;
  for (int component = 0; component < kDimension; ++component) {
    const double difference = first[component] - second[component];
    result += difference * difference;
  }
  return result;
}

double edge_length(const RegularMesh<kDimension>& mesh,
                   const GeometryBM& edge) {
  return std::sqrt(distance_squared(mesh.point(edge.vertex(0)),
                                    mesh.point(edge.vertex(1))));
}

double element_diameter(const RegularMesh<kDimension>& mesh,
                        const GeometryBM& element) {
  double diameter_squared = 0.0;
  for (int first = 0; first < element.n_vertex(); ++first) {
    for (int second = first + 1; second < element.n_vertex(); ++second) {
      diameter_squared = std::max(
          diameter_squared,
          distance_squared(mesh.point(element.vertex(first)),
                           mesh.point(element.vertex(second))));
    }
  }
  return std::sqrt(diameter_squared);
}

Point<kDimension> edge_midpoint(const RegularMesh<kDimension>& mesh,
                                const GeometryBM& edge) {
  Point<kDimension> midpoint;
  midpoint[0] = 0.5 * (mesh.point(edge.vertex(0))[0] +
                       mesh.point(edge.vertex(1))[0]);
  midpoint[1] = 0.5 * (mesh.point(edge.vertex(0))[1] +
                       mesh.point(edge.vertex(1))[1]);
  return midpoint;
}

void build_p1_space(
    RegularMesh<kDimension>& mesh,
    std::vector<TemplateElement<double, kDimension, kDimension> >&
        template_elements,
    FEMSpace<double, kDimension>& fem_space) {
  fem_space.reinit(mesh, template_elements);
  const int element_count = mesh.n_geometry(kDimension);
  fem_space.element().resize(element_count);
  for (int element_index = 0; element_index < element_count;
       ++element_index) {
    const int vertex_count =
        mesh.geometry(kDimension, element_index).n_vertex();
    if (vertex_count == 3) {
      fem_space.element(element_index).reinit(
          fem_space, element_index, 0);
    } else if (vertex_count == 4) {
      fem_space.element(element_index).reinit(
          fem_space, element_index, 1);
    } else {
      throw std::runtime_error(
          "Unsupported regular element: expected triangle or twin triangle.");
    }
  }
  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  // A corner vertex can carry the marker of either incident boundary in an
  // EasyMesh file.  The cooled bottom boundary has Dirichlet priority, so
  // assign its marker geometrically to every P1 dof on y=0.  In particular,
  // this prevents the bottom-right corner from being left unconstrained when
  // the mesh stores that vertex with the right-side Neumann marker.
  const double tolerance = 1.0e-12;
  for (unsigned int dof = 0; dof < fem_space.n_dof(); ++dof) {
    if (std::abs(fem_space.dofInfo(dof).interp_point[1]) < tolerance) {
      fem_space.dofBoundaryMark(dof) = kBottomBoundary;
    }
  }
}

Vector<double> assemble_flux_functional_gradient(
    RegularMesh<kDimension>& mesh,
    FEMSpace<double, kDimension>& fem_space) {
  Vector<double> gradient(fem_space.n_dof());
  const double tolerance = 1.0e-12;

  // On the bottom boundary n=(0,-1), hence the outward heat flux is
  // -k grad(u_h).n = k d_y u_h.  Since the P1 basis gradient is constant in
  // each triangle, midpoint integration is exact on every selected edge.
  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& geometry =
        mesh.geometry(kDimension, element_index);
    const Element<double, kDimension>& element =
        fem_space.element(element_index);
    for (int local_edge = 0; local_edge < geometry.n_boundary();
         ++local_edge) {
      const GeometryBM& edge = mesh.geometry(
          kDimension - 1, geometry.boundary(local_edge));
      if (edge.boundaryMark() != kBottomBoundary) {
        continue;
      }
      const Point<kDimension> midpoint = edge_midpoint(mesh, edge);
      if (midpoint[0] < kFluxPatchLeft - tolerance ||
          midpoint[0] > kFluxPatchRight + tolerance) {
        continue;
      }
      const double length = edge_length(mesh, edge);
      const std::vector<std::vector<double> > basis_gradients =
          element.basis_function_gradient(midpoint);
      for (unsigned int local_dof = 0;
           local_dof < element.dof().size(); ++local_dof) {
        const int global_dof = element.dof()[local_dof];
        gradient(global_dof) +=
            length * kConductivity * basis_gradients[local_dof][1];
      }
    }
  }
  return gradient;
}

struct ResidualFields {
  std::vector<double> eta_squared;
  std::vector<double> strong_residual_mean_squared;
  std::vector<double> cell_residual_squared;
  std::vector<double> interior_jump_squared;
  std::vector<double> neumann_boundary_squared;
  double estimator;
};

struct PrimalResult {
  ResidualFields residual;
  double flux_functional;
  double l2_error;
  int elements;
  int degrees_of_freedom;
};

ResidualFields compute_residual_fields(
    RegularMesh<kDimension>& mesh,
    FEMSpace<double, kDimension>& fem_space,
    const FEMFunction<double, kDimension>& solution) {
  const int element_count = mesh.n_geometry(kDimension);
  const int edge_count = mesh.n_geometry(kDimension - 1);
  ResidualFields result;
  result.eta_squared.assign(element_count, 0.0);
  result.strong_residual_mean_squared.assign(element_count, 0.0);
  result.cell_residual_squared.assign(element_count, 0.0);
  result.interior_jump_squared.assign(element_count, 0.0);
  result.neumann_boundary_squared.assign(element_count, 0.0);
  std::vector<std::vector<int> > edge_elements(edge_count);

  for (int element_index = 0; element_index < element_count;
       ++element_index) {
    const GeometryBM& geometry =
        mesh.geometry(kDimension, element_index);
    const Element<double, kDimension>& element =
        fem_space.element(element_index);
    const double diameter = element_diameter(mesh, geometry);
    const QuadratureInfo<kDimension>& quadrature =
        element.findQuadratureInfo(5);
    const std::vector<Point<kDimension> > points =
        element.local_to_global(quadrature.quadraturePoint());
    const std::vector<double> jacobians =
        element.local_to_global_jacobian(quadrature.quadraturePoint());
    const double reference_volume = element.templateElement().volume();
    double volume_residual_squared = 0.0;
    double volume = 0.0;
    for (int q = 0; q < quadrature.n_quadraturePoint(); ++q) {
      const double source = manufactured_heat_source(points[q]);
      const double weight = quadrature.weight(q) * jacobians[q] *
                            reference_volume;
      // For P1 elements, Delta T_h vanishes inside each triangle, hence the
      // elementwise strong residual Q + Delta T_h equals Q.
      volume_residual_squared += weight * source * source;
      volume += weight;
    }
    result.strong_residual_mean_squared[element_index] =
        volume > 0.0 ? volume_residual_squared / volume : 0.0;
    const double cell_contribution =
        diameter * diameter * volume_residual_squared;
    result.eta_squared[element_index] += cell_contribution;
    result.cell_residual_squared[element_index] = cell_contribution;

    for (int local_edge = 0; local_edge < geometry.n_boundary();
         ++local_edge) {
      edge_elements[geometry.boundary(local_edge)].push_back(element_index);
    }
  }

  for (int edge_index = 0; edge_index < edge_count; ++edge_index) {
    const GeometryBM& edge = mesh.geometry(kDimension - 1, edge_index);
    const double length = edge_length(mesh, edge);
    const Point<kDimension> midpoint = edge_midpoint(mesh, edge);
    const Point<kDimension>& first = mesh.point(edge.vertex(0));
    const Point<kDimension>& second = mesh.point(edge.vertex(1));
    const double normal_x = (second[1] - first[1]) / length;
    const double normal_y = -(second[0] - first[0]) / length;
    const std::vector<int>& neighbors = edge_elements[edge_index];
    if (neighbors.size() == 2) {
      const std::vector<double> gradient_left = solution.gradient(
          midpoint, fem_space.element(neighbors[0]));
      const std::vector<double> gradient_right = solution.gradient(
          midpoint, fem_space.element(neighbors[1]));
      const double jump = kConductivity *
          ((gradient_left[0] - gradient_right[0]) * normal_x +
           (gradient_left[1] - gradient_right[1]) * normal_y);
      const double contribution =
          0.5 * length * length * jump * jump;
      for (int neighbor : neighbors) {
        result.eta_squared[neighbor] += contribution;
        result.interior_jump_squared[neighbor] += contribution;
      }
    } else if (neighbors.size() == 1) {
      if (edge.boundaryMark() == kBottomBoundary) {
        continue;
      }
      const std::vector<double> gradient = solution.gradient(
          midpoint, fem_space.element(neighbors[0]));
      const double normal_flux = kConductivity *
          (gradient[0] * normal_x + gradient[1] * normal_y);
      const double contribution =
          length * length * normal_flux * normal_flux;
      result.eta_squared[neighbors[0]] += contribution;
      result.neumann_boundary_squared[neighbors[0]] += contribution;
    }
  }

  result.estimator = std::sqrt(std::accumulate(
      result.eta_squared.begin(), result.eta_squared.end(), 0.0));
  return result;
}

double evaluate_flux_functional(
    RegularMesh<kDimension>& mesh,
    FEMSpace<double, kDimension>& fem_space,
    const FEMFunction<double, kDimension>& solution) {
  double functional = 0.0;
  const double tolerance = 1.0e-12;
  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& geometry =
        mesh.geometry(kDimension, element_index);
    const Element<double, kDimension>& element =
        fem_space.element(element_index);
    for (int local_edge = 0; local_edge < geometry.n_boundary();
         ++local_edge) {
      const GeometryBM& edge = mesh.geometry(
          kDimension - 1, geometry.boundary(local_edge));
      if (edge.boundaryMark() != kBottomBoundary) {
        continue;
      }
      const Point<kDimension> midpoint = edge_midpoint(mesh, edge);
      if (midpoint[0] < kFluxPatchLeft - tolerance ||
          midpoint[0] > kFluxPatchRight + tolerance) {
        continue;
      }
      const std::vector<double> gradient =
          solution.gradient(midpoint, element);
      functional += edge_length(mesh, edge) *
                    kConductivity * gradient[1];
    }
  }
  return functional;
}

double compute_l2_error(
    FEMSpace<double, kDimension>& fem_space,
    const FEMFunction<double, kDimension>& solution) {
  double error_squared = 0.0;
  for (int element_index = 0;
       element_index < fem_space.mesh().n_geometry(kDimension);
       ++element_index) {
    const Element<double, kDimension>& element =
        fem_space.element(element_index);
    const QuadratureInfo<kDimension>& quadrature =
        element.findQuadratureInfo(5);
    const std::vector<Point<kDimension> > points =
        element.local_to_global(quadrature.quadraturePoint());
    const std::vector<double> jacobians =
        element.local_to_global_jacobian(quadrature.quadraturePoint());
    const std::vector<double> values = solution.value(points, element);
    const double reference_volume = element.templateElement().volume();
    for (int q = 0; q < quadrature.n_quadraturePoint(); ++q) {
      const double exact = exact_temperature(points[q]);
      const double difference = values[q] - exact;
      const double weight = quadrature.weight(q) * jacobians[q] *
                            reference_volume;
      error_squared += weight * difference * difference;
    }
  }
  return std::sqrt(error_squared);
}

PrimalResult solve_manufactured_primal(
    RegularMesh<kDimension>& mesh,
    std::vector<TemplateElement<double, kDimension, kDimension> >&
        template_elements) {
  FEMSpace<double, kDimension> fem_space;
  build_p1_space(mesh, template_elements, fem_space);

  StiffMatrix<kDimension, double> stiffness(fem_space);
  stiffness.algebricAccuracy() = 4;
  stiffness.build();
  Vector<double> right_hand_side;
  Operator::L2Discretize(
      &manufactured_heat_source, fem_space, right_hand_side, 5);
  FEMFunction<double, kDimension> solution(fem_space);
  BoundaryFunction<double, kDimension> cold_bottom(
      BoundaryConditionInfo::DIRICHLET,
      kBottomBoundary,
      &zero_boundary_value);
  BoundaryConditionAdmin<double, kDimension> boundary_admin(fem_space);
  boundary_admin.add(cold_bottom);
  boundary_admin.apply(stiffness, solution, right_hand_side);

  AMGSolver solver(stiffness);
  solver.solve(solution, right_hand_side, 1.0e-11, 500);
  solution.writeOpenDXData("manufactured_temperature.dx");

  PrimalResult result;
  result.residual = compute_residual_fields(mesh, fem_space, solution);
  result.flux_functional =
      evaluate_flux_functional(mesh, fem_space, solution);
  result.l2_error = compute_l2_error(fem_space, solution);
  result.elements = mesh.n_geometry(kDimension);
  result.degrees_of_freedom = fem_space.n_dof();
  return result;
}

struct DualResult {
  std::vector<double> element_rms_magnitude;
  double minimum;
  double maximum;
  int elements;
  int degrees_of_freedom;
};

void write_dual_nodal_data(
    FEMSpace<double, kDimension>& fem_space,
    const FEMFunction<double, kDimension>& dual,
    const std::string& filename) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  output << "# dof x y psi_h\n"
         << std::scientific << std::setprecision(12);
  for (unsigned int dof = 0; dof < fem_space.n_dof(); ++dof) {
    const Point<kDimension>& point = fem_space.dofInfo(dof).interp_point;
    output << dof << ' ' << point[0] << ' ' << point[1] << ' '
           << dual(dof) << '\n';
  }
}

DualResult solve_boundary_flux_dual(
    RegularMesh<kDimension>& mesh,
    std::vector<TemplateElement<double, kDimension, kDimension> >&
        template_elements) {
  FEMSpace<double, kDimension> fem_space;
  build_p1_space(mesh, template_elements, fem_space);

  StiffMatrix<kDimension, double> jacobian(fem_space);
  jacobian.algebricAccuracy() = 4;
  jacobian.build();

  // J_h(U) = integral_patch -k grad(u_h).n ds is linear in U.
  // The fully discrete dual is A_h^T psi_h = -grad_U J_h.
  Vector<double> dual_right_hand_side =
      assemble_flux_functional_gradient(mesh, fem_space);
  for (unsigned int dof = 0; dof < dual_right_hand_side.size(); ++dof) {
    dual_right_hand_side(dof) *= -1.0;
  }

  FEMFunction<double, kDimension> dual(fem_space);
  BoundaryFunction<double, kDimension> homogeneous_bottom(
      BoundaryConditionInfo::DIRICHLET,
      kBottomBoundary,
      &zero_boundary_value);
  BoundaryConditionAdmin<double, kDimension> boundary_admin(fem_space);
  boundary_admin.add(homogeneous_bottom);
  boundary_admin.apply(jacobian, dual, dual_right_hand_side);

  AMGSolver solver(jacobian);
  solver.solve(dual, dual_right_hand_side, 1.0e-11, 500);
  dual.writeOpenDXData("boundary_flux_dual_signed.dx");
  write_dual_nodal_data(
      fem_space, dual, "boundary_flux_dual_signed.dat");

  DualResult result;
  result.element_rms_magnitude.assign(
      mesh.n_geometry(kDimension), 0.0);
  result.minimum = 0.0;
  result.maximum = 0.0;
  result.elements = mesh.n_geometry(kDimension);
  result.degrees_of_freedom = fem_space.n_dof();
  for (unsigned int dof = 0; dof < dual.size(); ++dof) {
    result.minimum = std::min(result.minimum, dual(dof));
    result.maximum = std::max(result.maximum, dual(dof));
  }

  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const Element<double, kDimension>& element =
        fem_space.element(element_index);
    const QuadratureInfo<kDimension>& quadrature =
        element.findQuadratureInfo(4);
    const std::vector<Point<kDimension> > quadrature_points =
        element.local_to_global(quadrature.quadraturePoint());
    const std::vector<double> jacobians =
        element.local_to_global_jacobian(quadrature.quadraturePoint());
    const std::vector<double> values =
        dual.value(quadrature_points, element);
    const double reference_volume = element.templateElement().volume();
    double integral = 0.0;
    double volume = 0.0;
    for (int q = 0; q < quadrature.n_quadraturePoint(); ++q) {
      const double weight = quadrature.weight(q) * jacobians[q] *
                            reference_volume;
      integral += weight * values[q] * values[q];
      volume += weight;
    }
    result.element_rms_magnitude[element_index] =
        volume > 0.0 ? std::sqrt(integral / volume) : 0.0;
  }
  return result;
}

struct Bounds {
  double xmin;
  double xmax;
  double ymin;
  double ymax;
};

Bounds mesh_bounds(const RegularMesh<kDimension>& mesh) {
  Bounds bounds{mesh.point(0)[0], mesh.point(0)[0],
                mesh.point(0)[1], mesh.point(0)[1]};
  for (unsigned int point = 1; point < mesh.n_point(); ++point) {
    bounds.xmin = std::min(bounds.xmin, mesh.point(point)[0]);
    bounds.xmax = std::max(bounds.xmax, mesh.point(point)[0]);
    bounds.ymin = std::min(bounds.ymin, mesh.point(point)[1]);
    bounds.ymax = std::max(bounds.ymax, mesh.point(point)[1]);
  }
  return bounds;
}

double svg_x(double x, const Bounds& bounds) {
  return 55.0 + 690.0 * (x - bounds.xmin) /
                    (bounds.xmax - bounds.xmin);
}

double svg_y(double y, const Bounds& bounds) {
  return 745.0 - 690.0 * (y - bounds.ymin) /
                     (bounds.ymax - bounds.ymin);
}

void write_residual_svg(
    const RegularMesh<kDimension>& mesh,
    const std::vector<double>& squared_values,
    const std::string& filename,
    const std::string& title,
    double display_exponent) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  const Bounds bounds = mesh_bounds(mesh);
  const double maximum_squared =
      *std::max_element(squared_values.begin(), squared_values.end());
  const double maximum = std::sqrt(maximum_squared);
  output << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 800 820\">\n"
         << "<rect width=\"800\" height=\"820\" fill=\"white\"/>\n"
         << "<text x=\"400\" y=\"30\" text-anchor=\"middle\" "
         << "font-family=\"sans-serif\" font-size=\"20\">"
         << title << "</text>\n";
  output << std::setprecision(10);
  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& element =
        mesh.geometry(kDimension, element_index);
    const double magnitude = std::sqrt(squared_values[element_index]);
    const double normalized = maximum > 0.0 ? magnitude / maximum : 0.0;
    const double focused = std::pow(normalized, display_exponent);
    const int red = 250;
    const int green = static_cast<int>(246.0 - 155.0 * focused);
    const int blue = static_cast<int>(232.0 - 185.0 * focused);
    output << "<polygon points=\"";
    for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
      const Point<kDimension>& point =
          mesh.point(element.vertex(vertex));
      output << svg_x(point[0], bounds) << ','
             << svg_y(point[1], bounds) << ' ';
    }
    output << "\" fill=\"rgb(" << red << ',' << green << ',' << blue
           << ")\" stroke=\"#738496\" stroke-width=\"0.18\"/>\n";
  }

  const double patch_left = svg_x(kFluxPatchLeft, bounds);
  const double patch_right = svg_x(kFluxPatchRight, bounds);
  const double bottom = svg_y(0.0, bounds);
  output << "<line x1=\"" << patch_left << "\" y1=\"" << bottom
         << "\" x2=\"" << patch_right << "\" y2=\"" << bottom
         << "\" stroke=\"#0f766e\" stroke-width=\"7\"/>\n"
         << "<text x=\"" << 0.5 * (patch_left + patch_right)
         << "\" y=\"" << bottom - 12.0
         << "\" text-anchor=\"middle\" font-family=\"sans-serif\" "
         << "font-size=\"14\">target flux patch</text>\n"
         << "<text x=\"55\" y=\"785\" font-family=\"sans-serif\" "
         << "font-size=\"13\">max displayed magnitude = "
         << maximum << "</text>\n"
         << "<text x=\"745\" y=\"806\" text-anchor=\"end\" "
         << "font-family=\"sans-serif\" font-size=\"12\">"
         << "color exponent = " << display_exponent << "</text>\n"
         << "</svg>\n";
}

void write_dual_magnitude_svg(
    const RegularMesh<kDimension>& mesh,
    const DualResult& dual,
    const std::string& filename,
    double display_exponent,
    const std::string& scale_label) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  const Bounds bounds = mesh_bounds(mesh);
  const double maximum = *std::max_element(
      dual.element_rms_magnitude.begin(),
      dual.element_rms_magnitude.end());

  output << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 800 820\">\n"
         << "<rect width=\"800\" height=\"820\" fill=\"white\"/>\n"
         << "<text x=\"400\" y=\"30\" text-anchor=\"middle\" "
         << "font-family=\"sans-serif\" font-size=\"20\">"
         << "Discrete dual for heat flux through a cooled patch ("
         << scale_label << ")</text>\n";
  output << std::setprecision(10);

  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& element =
        mesh.geometry(kDimension, element_index);
    const double normalized = maximum > 0.0
        ? dual.element_rms_magnitude[element_index] / maximum
        : 0.0;
    // The exponent changes only display contrast, not the dual values.
    const double focused = std::pow(normalized, display_exponent);
    const int red = static_cast<int>(244.0 - 180.0 * focused);
    const int green = static_cast<int>(248.0 - 105.0 * focused);
    const int blue = static_cast<int>(246.0 - 65.0 * focused);
    output << "<polygon points=\"";
    for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
      const Point<kDimension>& point =
          mesh.point(element.vertex(vertex));
      output << svg_x(point[0], bounds) << ','
             << svg_y(point[1], bounds) << ' ';
    }
    output << "\" fill=\"rgb(" << red << ',' << green << ',' << blue
           << ")\" stroke=\"#738496\" stroke-width=\"0.18\"/>\n";
  }

  const double patch_left = svg_x(kFluxPatchLeft, bounds);
  const double patch_right = svg_x(kFluxPatchRight, bounds);
  const double bottom = svg_y(0.0, bounds);
  output << "<line x1=\"" << patch_left << "\" y1=\"" << bottom
         << "\" x2=\"" << patch_right << "\" y2=\"" << bottom
         << "\" stroke=\"#0f766e\" stroke-width=\"7\"/>\n"
         << "<line x1=\"" << patch_left << "\" y1=\"" << bottom - 9.0
         << "\" x2=\"" << patch_left << "\" y2=\"" << bottom + 3.0
         << "\" stroke=\"#0f766e\" stroke-width=\"2\"/>\n"
         << "<line x1=\"" << patch_right << "\" y1=\"" << bottom - 9.0
         << "\" x2=\"" << patch_right << "\" y2=\"" << bottom + 3.0
         << "\" stroke=\"#0f766e\" stroke-width=\"2\"/>\n"
         << "<text x=\"" << 0.5 * (patch_left + patch_right)
         << "\" y=\"" << bottom - 12.0
         << "\" text-anchor=\"middle\" font-family=\"sans-serif\" "
         << "font-size=\"14\">target flux patch</text>\n"
         << "<text x=\"55\" y=\"785\" font-family=\"sans-serif\" "
         << "font-size=\"13\">min psi_h = " << dual.minimum
         << ", max psi_h = " << dual.maximum << "</text>\n"
         << "<text x=\"745\" y=\"806\" text-anchor=\"end\" "
         << "font-family=\"sans-serif\" font-size=\"12\">"
         << "Dual location is determined by J, not by Q; color exponent = "
         << display_exponent << ".</text>\n"
         << "</svg>\n";
}

void write_dual_magnitude_data(
    const RegularMesh<kDimension>& mesh,
    const DualResult& dual,
    const std::string& filename) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  output << "# element centroid_x centroid_y dual_rms_magnitude\n"
         << std::scientific << std::setprecision(12);
  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& element =
        mesh.geometry(kDimension, element_index);
    double x = 0.0;
    double y = 0.0;
    for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
      x += mesh.point(element.vertex(vertex))[0];
      y += mesh.point(element.vertex(vertex))[1];
    }
    x /= element.n_vertex();
    y /= element.n_vertex();
    output << element_index << ' ' << x << ' ' << y << ' '
           << dual.element_rms_magnitude[element_index] << '\n';
  }
}

void write_residual_data(
    const RegularMesh<kDimension>& mesh,
    const ResidualFields& residual,
    const std::string& filename) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  output << "# element centroid_x centroid_y eta_squared strong_rms "
         << "cell_squared jump_squared neumann_squared\n"
         << std::scientific << std::setprecision(12);
  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& element =
        mesh.geometry(kDimension, element_index);
    double x = 0.0;
    double y = 0.0;
    for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
      x += mesh.point(element.vertex(vertex))[0];
      y += mesh.point(element.vertex(vertex))[1];
    }
    x /= element.n_vertex();
    y /= element.n_vertex();
    output << element_index << ' ' << x << ' ' << y << ' '
           << residual.eta_squared[element_index] << ' '
           << std::sqrt(
                  residual.strong_residual_mean_squared[element_index])
           << ' ' << residual.cell_residual_squared[element_index] << ' '
           << residual.interior_jump_squared[element_index] << ' '
           << residual.neumann_boundary_squared[element_index] << '\n';
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  const std::string root_mesh = argc >= 2 ? argv[1] : "D";
  try {
    TemplateGeometry<kDimension> triangle_geometry;
    triangle_geometry.readData("triangle.tmp_geo");
    CoordTransform<kDimension, kDimension> triangle_transform;
    triangle_transform.readData("triangle.crd_trs");
    TemplateDOF<kDimension> triangle_dof(triangle_geometry);
    triangle_dof.readData("triangle.1.tmp_dof");
    BasisFunctionAdmin<double, kDimension, kDimension> triangle_basis(
        triangle_dof);
    triangle_basis.readData("triangle.1.bas_fun");

    TemplateGeometry<kDimension> twin_geometry;
    twin_geometry.readData("twin_triangle.tmp_geo");
    CoordTransform<kDimension, kDimension> twin_transform;
    twin_transform.readData("twin_triangle.crd_trs");
    TemplateDOF<kDimension> twin_dof(twin_geometry);
    twin_dof.readData("twin_triangle.1.tmp_dof");
    BasisFunctionAdmin<double, kDimension, kDimension> twin_basis(twin_dof);
    twin_basis.readData("twin_triangle.1.bas_fun");

    std::vector<TemplateElement<double, kDimension, kDimension> >
        template_elements(2);
    template_elements[0].reinit(triangle_geometry,
                                triangle_dof,
                                triangle_transform,
                                triangle_basis);
    template_elements[1].reinit(twin_geometry,
                                twin_dof,
                                twin_transform,
                                twin_basis);

    HGeometryTree<kDimension> hierarchy;
    hierarchy.readEasyMesh(root_mesh);
    IrregularMesh<kDimension> irregular_mesh(hierarchy);
    irregular_mesh.globalRefine(kUniformRefinementLevel);
    irregular_mesh.semiregularize();
    irregular_mesh.regularize(false);
    RegularMesh<kDimension>& mesh = irregular_mesh.regularMesh();

    const PrimalResult primal =
        solve_manufactured_primal(mesh, template_elements);
    const DualResult dual =
        solve_boundary_flux_dual(mesh, template_elements);
    write_residual_svg(mesh,
                       primal.residual.strong_residual_mean_squared,
                       "manufactured_strong_residual.svg",
                       "Strong residual RMS: Q + Delta T_h = Q",
                       1.0);
    write_residual_svg(mesh,
                       primal.residual.eta_squared,
                       "manufactured_full_residual_indicator.svg",
                       "Full residual indicator: cell plus flux jumps",
                       1.0);
    write_residual_data(
        mesh, primal.residual, "manufactured_residual.dat");
    write_dual_magnitude_svg(mesh,
                             dual,
                             "boundary_flux_dual_magnitude.svg",
                             1.0,
                             "linear color scale");
    write_dual_magnitude_svg(mesh,
                             dual,
                             "boundary_flux_dual_magnitude_focused.svg",
                             2.0,
                             "focused color scale");
    write_dual_magnitude_data(
        mesh, dual, "boundary_flux_dual_magnitude.dat");

    std::cout << std::setprecision(8)
              << "Boundary-flux discrete dual\n"
              << "  target patch        : [" << kFluxPatchLeft << ", "
              << kFluxPatchRight << "] x {0}\n"
              << "  manufactured alpha  : "
              << kManufacturedAmplitude << '\n'
              << "  uniform mesh level  : " << kUniformRefinementLevel
              << '\n'
              << "  elements            : " << dual.elements << '\n'
              << "  dofs                : " << dual.degrees_of_freedom << '\n'
              << "  global estimator    : "
              << primal.residual.estimator << '\n'
              << "  L2 temperature error: " << primal.l2_error << '\n'
              << "  numerical J_h       : "
              << primal.flux_functional << '\n'
              << "  exact J             : "
              << exact_flux_functional() << '\n'
              << "  |J_h-J|             : "
              << std::abs(primal.flux_functional -
                          exact_flux_functional()) << '\n'
              << "  dual range          : [" << dual.minimum << ", "
              << dual.maximum << "]\n"
              << "  note                : moving the source does not change "
              << "this dual\n";
  } catch (const std::exception& error) {
    std::cerr << "boundary_flux_dual: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
