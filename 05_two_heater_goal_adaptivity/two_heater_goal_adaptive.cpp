#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
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

namespace {

constexpr int kDimension = 2;
constexpr double kConductivity = 1.0;
constexpr int kBoundaryMarks[4] = {1, 2, 3, 4};
constexpr double kDistractorStrength = 240.0;
constexpr double kDistractorCenterX = 0.15;
constexpr double kDistractorCenterY = 0.85;
constexpr double kDistractorSigmaX = 0.050;
constexpr double kDistractorSigmaY = 0.050;
constexpr double kTargetStrength = 25.0;
constexpr double kTargetCenterX = 0.55;
constexpr double kTargetCenterY = 0.40;
constexpr double kTargetSigma = 0.080;
constexpr double kSensorCenterX = 0.82;
constexpr double kSensorCenterY = 0.18;
constexpr double kSensorSigma = 0.120;
constexpr int kDefaultAdaptationRounds = 6;
constexpr double kRefineFraction = 0.05;
constexpr int kReferenceModes = 160;
constexpr int kReferenceIntervals = 4096;

void prepare_output_directories() {
  for (const char* directory : {
           "summary",
           "fields/problem",
           "fields/residual",
           "fields/dual",
           "fields/dwr",
           "meshes/residual",
           "meshes/dual",
           "meshes/dwr"}) {
    std::filesystem::create_directories(directory);
  }
}

double heat_source(const double* point) {
  const double distractor_dx = point[0] - kDistractorCenterX;
  const double distractor_dy = point[1] - kDistractorCenterY;
  const double target_dx = point[0] - kTargetCenterX;
  const double target_dy = point[1] - kTargetCenterY;
  const double distractor = kDistractorStrength * std::exp(
      -(distractor_dx * distractor_dx /
            (2.0 * kDistractorSigmaX * kDistractorSigmaX) +
        distractor_dy * distractor_dy /
            (2.0 * kDistractorSigmaY * kDistractorSigmaY)));
  const double target = kTargetStrength * std::exp(
      -(target_dx * target_dx + target_dy * target_dy) /
      (2.0 * kTargetSigma * kTargetSigma));
  return distractor + target;
}

double heat_source_at_point(const Point<kDimension>& point) {
  const double coordinates[kDimension] = {point[0], point[1]};
  return heat_source(coordinates);
}

double cold_temperature(const double*) {
  return 0.0;
}

double gaussian_interval_integral(double center, double sigma) {
  const double sqrt_two = std::sqrt(2.0);
  const double sqrt_pi_over_two = std::sqrt(std::acos(-1.0) / 2.0);
  return sigma * sqrt_pi_over_two *
         (std::erf((1.0 - center) / (sqrt_two * sigma)) -
          std::erf(-center / (sqrt_two * sigma)));
}

double sensor_normalization() {
  return gaussian_interval_integral(kSensorCenterX, kSensorSigma) *
         gaussian_interval_integral(kSensorCenterY, kSensorSigma);
}

double sensor_weight(const Point<kDimension>& point) {
  const double dx = point[0] - kSensorCenterX;
  const double dy = point[1] - kSensorCenterY;
  return std::exp(-(dx * dx + dy * dy) /
                  (2.0 * kSensorSigma * kSensorSigma)) /
         sensor_normalization();
}

double sensor_functional_gradient(const double* point) {
  const double dx = point[0] - kSensorCenterX;
  const double dy = point[1] - kSensorCenterY;
  return std::exp(-(dx * dx + dy * dy) /
                  (2.0 * kSensorSigma * kSensorSigma)) /
         sensor_normalization();
}

double gaussian_basis_integral(double center,
                               double sigma,
                               int mode,
                               int intervals) {
  if (intervals <= 0 || intervals % 2 != 0) {
    throw std::invalid_argument(
        "Composite Simpson quadrature needs a positive even interval count.");
  }
  const double pi = std::acos(-1.0);
  const double step = 1.0 / intervals;
  double sum = 0.0;
  for (int i = 0; i <= intervals; ++i) {
    const double coordinate = i * step;
    const double gaussian = std::exp(
        -(coordinate - center) * (coordinate - center) /
        (2.0 * sigma * sigma));
    const double basis = std::sqrt(2.0) *
                         std::sin(mode * pi * coordinate);
    const double simpson_weight =
        (i == 0 || i == intervals) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
    sum += simpson_weight * gaussian * basis;
  }
  return step * sum / 3.0;
}

double spectral_reference_sensor_temperature(int modes, int intervals) {
  const double pi = std::acos(-1.0);
  const double sensor_x_normalization =
      gaussian_interval_integral(kSensorCenterX, kSensorSigma);
  const double sensor_y_normalization =
      gaussian_interval_integral(kSensorCenterY, kSensorSigma);
  std::vector<double> source_x(modes, 0.0);
  std::vector<double> source_y(modes, 0.0);
  std::vector<double> target_source_x(modes, 0.0);
  std::vector<double> target_source_y(modes, 0.0);
  std::vector<double> sensor_x(modes, 0.0);
  std::vector<double> sensor_y(modes, 0.0);

  for (int mode = 0; mode < modes; ++mode) {
    const int mode_number = mode + 1;
    source_x[mode] = gaussian_basis_integral(
        kDistractorCenterX, kDistractorSigmaX, mode_number, intervals);
    source_y[mode] = gaussian_basis_integral(
        kDistractorCenterY, kDistractorSigmaY, mode_number, intervals);
    target_source_x[mode] = gaussian_basis_integral(
        kTargetCenterX, kTargetSigma, mode_number, intervals);
    target_source_y[mode] = gaussian_basis_integral(
        kTargetCenterY, kTargetSigma, mode_number, intervals);
    sensor_x[mode] = gaussian_basis_integral(
        kSensorCenterX, kSensorSigma, mode_number, intervals) /
        sensor_x_normalization;
    sensor_y[mode] = gaussian_basis_integral(
        kSensorCenterY, kSensorSigma, mode_number, intervals) /
        sensor_y_normalization;
  }

  double functional = 0.0;
  for (int x_mode = 0; x_mode < modes; ++x_mode) {
    for (int y_mode = 0; y_mode < modes; ++y_mode) {
      const double x_wave_number = (x_mode + 1) * pi;
      const double y_wave_number = (y_mode + 1) * pi;
      const double eigenvalue =
          x_wave_number * x_wave_number + y_wave_number * y_wave_number;
      const double source_coefficient =
          kDistractorStrength * source_x[x_mode] * source_y[y_mode] +
          kTargetStrength * target_source_x[x_mode] *
              target_source_y[y_mode];
      functional += source_coefficient *
          sensor_x[x_mode] * sensor_y[y_mode] / eigenvalue;
    }
  }
  return functional;
}

double distance_squared(const Point<kDimension>& first,
                        const Point<kDimension>& second) {
  const double dx = first[0] - second[0];
  const double dy = first[1] - second[1];
  return dx * dx + dy * dy;
}

double element_diameter(const RegularMesh<kDimension>& mesh,
                        const GeometryBM& geometry) {
  double diameter_squared = 0.0;
  for (int i = 0; i < geometry.n_vertex(); ++i) {
    for (int j = i + 1; j < geometry.n_vertex(); ++j) {
      diameter_squared = std::max(
          diameter_squared,
          distance_squared(mesh.point(geometry.vertex(i)),
                           mesh.point(geometry.vertex(j))));
    }
  }
  return std::sqrt(diameter_squared);
}

double edge_length(const RegularMesh<kDimension>& mesh,
                   const GeometryBM& edge) {
  if (edge.n_vertex() != 2) {
    throw std::runtime_error("Expected a two-vertex edge.");
  }
  return std::sqrt(distance_squared(mesh.point(edge.vertex(0)),
                                    mesh.point(edge.vertex(1))));
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

Point<kDimension> element_centroid(
    const RegularMesh<kDimension>& mesh,
    const GeometryBM& element) {
  Point<kDimension> centroid;
  centroid[0] = 0.0;
  centroid[1] = 0.0;
  for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
    centroid[0] += mesh.point(element.vertex(vertex))[0];
    centroid[1] += mesh.point(element.vertex(vertex))[1];
  }
  centroid[0] /= element.n_vertex();
  centroid[1] /= element.n_vertex();
  return centroid;
}

std::vector<double> unit_normal(const RegularMesh<kDimension>& mesh,
                                const GeometryBM& edge) {
  const Point<kDimension>& first = mesh.point(edge.vertex(0));
  const Point<kDimension>& second = mesh.point(edge.vertex(1));
  const double length = edge_length(mesh, edge);
  std::vector<double> normal(2);
  normal[0] = (second[1] - first[1]) / length;
  normal[1] = -(second[0] - first[0]) / length;
  return normal;
}

std::vector<double> outward_unit_normal(
    const RegularMesh<kDimension>& mesh,
    const GeometryBM& edge,
    const GeometryBM& element) {
  std::vector<double> normal = unit_normal(mesh, edge);
  const Point<kDimension> midpoint = edge_midpoint(mesh, edge);
  const Point<kDimension> centroid = element_centroid(mesh, element);
  const double toward_element =
      normal[0] * (centroid[0] - midpoint[0]) +
      normal[1] * (centroid[1] - midpoint[1]);
  if (toward_element > 0.0) {
    normal[0] *= -1.0;
    normal[1] *= -1.0;
  }
  return normal;
}

void build_p1_space(
    RegularMesh<kDimension>& mesh,
    std::vector<TemplateElement<double, kDimension, kDimension> >&
        template_elements,
    FEMSpace<double, kDimension>& fem_space) {
  fem_space.reinit(mesh, template_elements);
  const int element_count = mesh.n_geometry(kDimension);
  fem_space.element().resize(element_count);

  for (int i = 0; i < element_count; ++i) {
    const int vertex_count = mesh.geometry(kDimension, i).n_vertex();
    if (vertex_count == 3) {
      fem_space.element(i).reinit(fem_space, i, 0);
    } else if (vertex_count == 4) {
      fem_space.element(i).reinit(fem_space, i, 1);
    } else {
      throw std::runtime_error(
          "Unsupported regular element: expected triangle or twin triangle.");
    }
  }

  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();
}

void apply_zero_dirichlet(
    FEMSpace<double, kDimension>& fem_space,
    StiffMatrix<kDimension, double>& stiffness,
    FEMFunction<double, kDimension>& solution,
    Vector<double>& right_hand_side) {
  BoundaryFunction<double, kDimension> boundary_1(
      BoundaryConditionInfo::DIRICHLET,
      kBoundaryMarks[0],
      &cold_temperature);
  BoundaryFunction<double, kDimension> boundary_2(
      BoundaryConditionInfo::DIRICHLET,
      kBoundaryMarks[1],
      &cold_temperature);
  BoundaryFunction<double, kDimension> boundary_3(
      BoundaryConditionInfo::DIRICHLET,
      kBoundaryMarks[2],
      &cold_temperature);
  BoundaryFunction<double, kDimension> boundary_4(
      BoundaryConditionInfo::DIRICHLET,
      kBoundaryMarks[3],
      &cold_temperature);
  BoundaryConditionAdmin<double, kDimension> boundary_admin(fem_space);
  boundary_admin.add(boundary_1);
  boundary_admin.add(boundary_2);
  boundary_admin.add(boundary_3);
  boundary_admin.add(boundary_4);
  boundary_admin.apply(stiffness, solution, right_hand_side);
}

std::vector<double> compute_residual_indicator(
    RegularMesh<kDimension>& mesh,
    FEMSpace<double, kDimension>& fem_space,
    const FEMFunction<double, kDimension>& solution,
    std::vector<double>& strong_residual_mean_squared,
    std::vector<double>& cell_residual_squared,
    std::vector<double>& interior_jump_squared,
    std::vector<double>& neumann_boundary_squared) {
  const int element_count = mesh.n_geometry(kDimension);
  const int edge_count = mesh.n_geometry(kDimension - 1);
  std::vector<double> eta_squared(element_count, 0.0);
  strong_residual_mean_squared.assign(element_count, 0.0);
  cell_residual_squared.assign(element_count, 0.0);
  interior_jump_squared.assign(element_count, 0.0);
  neumann_boundary_squared.assign(element_count, 0.0);
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
    const std::vector<Point<kDimension> > quadrature_points =
        element.local_to_global(quadrature.quadraturePoint());
    const std::vector<double> jacobians =
        element.local_to_global_jacobian(quadrature.quadraturePoint());
    const double reference_volume = element.templateElement().volume();

    double volume_residual_squared = 0.0;
    double element_volume = 0.0;
    for (int q = 0; q < quadrature.n_quadraturePoint(); ++q) {
      const double source = heat_source(quadrature_points[q]);
      const double weight = quadrature.weight(q) * jacobians[q] *
                            reference_volume;
      // For conforming P1 elements and elementwise constant conductivity,
      // div(kappa grad T_h) vanishes inside each simplex.  The remaining
      // strong element residual is therefore the heat source.
      volume_residual_squared += weight * source * source;
      element_volume += weight;
    }
    strong_residual_mean_squared[element_index] =
        element_volume > 0.0 ? volume_residual_squared / element_volume : 0.0;
    const double cell_contribution =
        diameter * diameter * volume_residual_squared;
    eta_squared[element_index] += cell_contribution;
    cell_residual_squared[element_index] += cell_contribution;

    // A twin triangle is the union of template triangles (0,1,2) and
    // (0,2,3).  Their common edge (0,2) is internal to the macro element,
    // so it does not appear in RegularMesh::geometry(1).  Add its flux-jump
    // contribution explicitly.
    if (geometry.n_vertex() == 4) {
      const Point<kDimension>& vertex_0 =
          mesh.point(geometry.vertex(0));
      const Point<kDimension>& vertex_1 =
          mesh.point(geometry.vertex(1));
      const Point<kDimension>& vertex_2 =
          mesh.point(geometry.vertex(2));
      const Point<kDimension>& vertex_3 =
          mesh.point(geometry.vertex(3));
      Point<kDimension> centroid_upper;
      Point<kDimension> centroid_lower;
      for (int component = 0; component < kDimension; ++component) {
        centroid_upper[component] =
            (vertex_0[component] + vertex_1[component] +
             vertex_2[component]) /
            3.0;
        centroid_lower[component] =
            (vertex_0[component] + vertex_2[component] +
             vertex_3[component]) /
            3.0;
      }
      const std::vector<double> gradient_upper =
          solution.gradient(centroid_upper, element);
      const std::vector<double> gradient_lower =
          solution.gradient(centroid_lower, element);
      const double internal_length =
          std::sqrt(distance_squared(vertex_0, vertex_2));
      const double normal_x =
          (vertex_2[1] - vertex_0[1]) / internal_length;
      const double normal_y =
          -(vertex_2[0] - vertex_0[0]) / internal_length;
      const double internal_jump =
          kConductivity *
          ((gradient_upper[0] - gradient_lower[0]) * normal_x +
           (gradient_upper[1] - gradient_lower[1]) * normal_y);
      const double contribution =
          internal_length * internal_length *
          internal_jump * internal_jump;
      eta_squared[element_index] += contribution;
      interior_jump_squared[element_index] += contribution;
    }

    for (int local_edge = 0; local_edge < geometry.n_boundary();
         ++local_edge) {
      const int edge_index = geometry.boundary(local_edge);
      edge_elements[edge_index].push_back(element_index);
    }
  }

  for (int edge_index = 0; edge_index < edge_count; ++edge_index) {
    const GeometryBM& edge = mesh.geometry(kDimension - 1, edge_index);
    const double length = edge_length(mesh, edge);
    const Point<kDimension> midpoint = edge_midpoint(mesh, edge);
    const std::vector<double> normal = unit_normal(mesh, edge);
    const std::vector<int>& neighbors = edge_elements[edge_index];

    if (neighbors.size() == 2) {
      const std::vector<double> gradient_left =
          solution.gradient(midpoint, fem_space.element(neighbors[0]));
      const std::vector<double> gradient_right =
          solution.gradient(midpoint, fem_space.element(neighbors[1]));
      const double flux_jump =
          kConductivity *
          ((gradient_left[0] - gradient_right[0]) * normal[0] +
           (gradient_left[1] - gradient_right[1]) * normal[1]);
      const double contribution =
          0.5 * length * length * flux_jump * flux_jump;
      eta_squared[neighbors[0]] += contribution;
      eta_squared[neighbors[1]] += contribution;
      interior_jump_squared[neighbors[0]] += contribution;
      interior_jump_squared[neighbors[1]] += contribution;
    } else if (neighbors.size() == 1) {
      if (edge.boundaryMark() != 0) {
        continue;
      }
      const std::vector<double> gradient =
          solution.gradient(midpoint, fem_space.element(neighbors[0]));
      const double normal_flux =
          kConductivity *
          (gradient[0] * normal[0] + gradient[1] * normal[1]);
      const double contribution =
          length * length * normal_flux * normal_flux;
      eta_squared[neighbors[0]] += contribution;
      neumann_boundary_squared[neighbors[0]] += contribution;
    } else {
      throw std::runtime_error(
          "Unexpected edge adjacency while computing the residual.");
    }
  }

  return eta_squared;
}

struct ProportionalMarking {
  std::vector<bool> marked;
  double tolerance;
};

ProportionalMarking mark_largest_fraction(
    const std::vector<double>& indicator,
    double fraction) {
  if (indicator.empty() || fraction <= 0.0 || fraction > 1.0) {
    throw std::invalid_argument(
        "The refinement fraction must lie in (0,1].");
  }
  std::vector<int> order(indicator.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&indicator](int first, int second) {
              return indicator[first] > indicator[second];
            });

  const int selected = std::min(
      static_cast<int>(indicator.size()),
      std::max(1, static_cast<int>(
                      std::ceil(fraction * indicator.size()))));
  const double cutoff = indicator[order[selected - 1]];
  // MeshAdaptor refines on a strict comparison.  Move the quantile by one
  // floating-point value so the element at the cutoff is included.
  const double tolerance =
      std::nextafter(cutoff, -std::numeric_limits<double>::infinity());
  std::vector<bool> marked(indicator.size(), false);
  for (int i = 0; i < selected; ++i) {
    marked[order[i]] = true;
  }
  return ProportionalMarking{marked, tolerance};
}

void adapt_with_afepack_quantile(
    IrregularMesh<kDimension>& irregular_mesh,
    RegularMesh<kDimension>& regular_mesh,
    const std::vector<double>& values,
    double tolerance) {
  Indicator<kDimension> indicator(regular_mesh);
  std::copy(values.begin(), values.end(), indicator.begin());

  MeshAdaptor<kDimension> adaptor(irregular_mesh);
  adaptor.convergenceOrder() = 0.0;
  adaptor.refineStep() = 0;
  // In two dimensions the overflow test is
  // indicator > refineThreshold * 4 * tolerence.  The value 0.25 therefore
  // makes the AFEPack tolerance exactly the selected indicator quantile.
  adaptor.refineThreshold() = 0.25;
  adaptor.tolerence() = tolerance;
  adaptor.is_refine_only() = true;
  adaptor.setIndicator(indicator);
  adaptor.adapt();
}


std::vector<double> sample_at_element_centroids(
    const RegularMesh<kDimension>& mesh,
    double (*field)(const Point<kDimension>&)) {
  std::vector<double> values(mesh.n_geometry(kDimension), 0.0);
  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& element = mesh.geometry(kDimension, element_index);
    Point<kDimension> centroid;
    centroid[0] = 0.0;
    centroid[1] = 0.0;
    for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
      centroid[0] += mesh.point(element.vertex(vertex))[0];
      centroid[1] += mesh.point(element.vertex(vertex))[1];
    }
    centroid[0] /= element.n_vertex();
    centroid[1] /= element.n_vertex();
    values[element_index] = field(centroid);
  }
  return values;
}

void write_cell_field_data(const RegularMesh<kDimension>& mesh,
                           const std::vector<double>& values,
                           const std::string& filename) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  output << "# element centroid_x centroid_y value\n"
         << std::scientific << std::setprecision(12);
  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& element = mesh.geometry(kDimension, element_index);
    double x = 0.0;
    double y = 0.0;
    for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
      const Point<kDimension>& point = mesh.point(element.vertex(vertex));
      x += point[0];
      y += point[1];
    }
    output << element_index << ' '
           << x / element.n_vertex() << ' '
           << y / element.n_vertex() << ' '
           << values[element_index] << '\n';
  }
}

void write_indicator_data(const RegularMesh<kDimension>& mesh,
                          const std::vector<double>& eta_squared,
                          const std::vector<double>& strong_residual_squared,
                          const std::vector<double>& cell_residual_squared,
                          const std::vector<double>& interior_jump_squared,
                          const std::vector<double>& neumann_boundary_squared,
                          const std::vector<bool>& marked,
                          const std::string& filename) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  output << "# element centroid_x centroid_y eta_squared "
         << "strong_residual_mean_squared cell_squared "
         << "interior_jump_squared neumann_boundary_squared marked\n";
  output << std::scientific << std::setprecision(12);
  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& element = mesh.geometry(kDimension, element_index);
    double x = 0.0;
    double y = 0.0;
    for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
      const Point<kDimension>& point = mesh.point(element.vertex(vertex));
      x += point[0];
      y += point[1];
    }
    x /= element.n_vertex();
    y /= element.n_vertex();
    output << element_index << ' ' << x << ' ' << y << ' '
           << eta_squared[element_index] << ' '
           << strong_residual_squared[element_index] << ' '
           << cell_residual_squared[element_index] << ' '
           << interior_jump_squared[element_index] << ' '
           << neumann_boundary_squared[element_index] << ' '
           << (marked[element_index] ? 1 : 0) << '\n';
  }
}

void write_dual_magnitude_data(
    const RegularMesh<kDimension>& mesh,
    const std::vector<double>& magnitude,
    const std::vector<bool>& marked,
    const std::string& filename) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  output << "# element centroid_x centroid_y dual_rms marked\n"
         << std::scientific << std::setprecision(12);
  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& element = mesh.geometry(kDimension, element_index);
    double x = 0.0;
    double y = 0.0;
    for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
      const Point<kDimension>& point = mesh.point(element.vertex(vertex));
      x += point[0];
      y += point[1];
    }
    x /= element.n_vertex();
    y /= element.n_vertex();
    output << element_index << ' ' << x << ' ' << y << ' '
           << magnitude[element_index] << ' '
           << (marked[element_index] ? 1 : 0) << '\n';
  }
}

void write_discrete_dual_data(
    FEMSpace<double, kDimension>& fem_space,
    const FEMFunction<double, kDimension>& dual,
    const std::string& filename) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  output << "# dof x y psi_h\n" << std::scientific << std::setprecision(12);
  for (unsigned int i = 0; i < fem_space.n_dof(); ++i) {
    const Point<kDimension>& point = fem_space.dofInfo(i).interp_point;
    output << i << ' ' << point[0] << ' ' << point[1] << ' '
           << dual(i) << '\n';
  }
}

struct DiscreteDualResult {
  std::vector<double> element_rms_magnitude;
  double minimum;
  double maximum;
  int degrees_of_freedom;
  int elements;
};

struct SolveResult {
  std::vector<double> eta_squared;
  std::vector<double> strong_residual_mean_squared;
  std::vector<double> cell_residual_squared;
  std::vector<double> interior_jump_squared;
  std::vector<double> neumann_boundary_squared;
  double estimator;
  double sensor_temperature;
  int degrees_of_freedom;
};

struct AdaptationReport {
  int round;
  int elements_before;
  int dofs_before;
  double estimator;
  double sensor_temperature;
  int marked_elements;
  int elements_after;
  int twin_triangles_after;
};

struct DualAdaptationReport {
  int round;
  int elements_before;
  int dofs_before;
  double sensor_temperature;
  double minimum;
  double maximum;
  int marked_elements;
  int elements_after;
};

struct DWRResult {
  std::vector<double> signed_element_contribution;
  std::vector<double> absolute_element_indicator;
  std::vector<double> primal_residual_norm;
  std::vector<double> dual_solution_rms;
  std::vector<double> dual_correction_rms;
  std::vector<double> residual_dual_correction_product;
  std::vector<double> residual_absolute_dual_product;
  double signed_estimate;
  double absolute_sum;
  double sensor_temperature;
  int primal_degrees_of_freedom;
  int enriched_dual_degrees_of_freedom;
};

struct DWRAdaptationReport {
  int round;
  int elements_before;
  int dofs_before;
  double sensor_temperature;
  double signed_estimate;
  double absolute_sum;
  int enriched_dual_dofs;
  int marked_elements;
  int elements_after;
};

struct FunctionalRecord {
  int level;
  int elements;
  int degrees_of_freedom;
  double value;
};

int count_twin_triangles(const RegularMesh<kDimension>& mesh) {
  int count = 0;
  for (unsigned int i = 0; i < mesh.n_geometry(kDimension); ++i) {
    if (mesh.geometry(kDimension, i).n_vertex() == 4) {
      ++count;
    }
  }
  return count;
}

double evaluate_sensor_temperature(
    FEMSpace<double, kDimension>& fem_space,
    const FEMFunction<double, kDimension>& solution) {
  double sensor_temperature = 0.0;
  for (int element_index = 0;
       element_index < fem_space.mesh().n_geometry(kDimension);
       ++element_index) {
    const Element<double, kDimension>& element =
        fem_space.element(element_index);
    const QuadratureInfo<kDimension>& quadrature =
        element.findQuadratureInfo(5);
    const std::vector<Point<kDimension> > quadrature_points =
        element.local_to_global(quadrature.quadraturePoint());
    const std::vector<double> jacobians =
        element.local_to_global_jacobian(quadrature.quadraturePoint());
    const std::vector<double> temperatures =
        solution.value(quadrature_points, element);
    const double reference_volume = element.templateElement().volume();

    for (int q = 0; q < quadrature.n_quadraturePoint(); ++q) {
      const double weight = quadrature.weight(q) * jacobians[q] *
                            reference_volume;
      sensor_temperature += weight * sensor_weight(quadrature_points[q]) *
                            temperatures[q];
    }
  }
  return sensor_temperature;
}

DiscreteDualResult solve_discrete_dual(
    RegularMesh<kDimension>& mesh,
    std::vector<TemplateElement<double, kDimension, kDimension> >&
        template_elements,
    const std::string& output_stem) {
  FEMSpace<double, kDimension> fem_space;
  build_p1_space(mesh, template_elements, fem_space);

  // Assemble the Jacobian dR_h/dU_h of the discrete residual
  // R_h(U) = A U - b.  For this Poisson P1 discretization A is symmetric,
  // hence the same matrix represents A^T.  This is the only Poisson-specific
  // shortcut in the fully discrete adjoint construction.
  StiffMatrix<kDimension, double> discrete_jacobian(fem_space);
  discrete_jacobian.algebricAccuracy() = 4;
  discrete_jacobian.build();

  // Assemble g_h = grad_U J_h directly in the discrete basis, then solve
  // (dR_h/dU_h)^T psi_h = -g_h.
  Vector<double> dual_right_hand_side;
  Operator::L2Discretize(&sensor_functional_gradient,
                         fem_space,
                         dual_right_hand_side,
                         5);
  for (unsigned int i = 0; i < dual_right_hand_side.size(); ++i) {
    dual_right_hand_side(i) *= -1.0;
  }

  FEMFunction<double, kDimension> dual(fem_space);
  apply_zero_dirichlet(fem_space,
                       discrete_jacobian,
                       dual,
                       dual_right_hand_side);

  AMGSolver dual_solver(discrete_jacobian);
  dual_solver.solve(dual, dual_right_hand_side, 1.0e-11, 500);

  DiscreteDualResult result;
  result.element_rms_magnitude.assign(
      mesh.n_geometry(kDimension), 0.0);
  result.minimum = 0.0;
  result.maximum = 0.0;
  result.degrees_of_freedom = fem_space.n_dof();
  result.elements = mesh.n_geometry(kDimension);
  for (unsigned int i = 0; i < dual.size(); ++i) {
    result.minimum = std::min(result.minimum, dual(i));
    result.maximum = std::max(result.maximum, dual(i));
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

  dual.writeOpenDXData(output_stem + ".dx");
  write_discrete_dual_data(fem_space,
                           dual,
                           output_stem + ".dat");

  return result;
}

SolveResult solve_temperature(
    RegularMesh<kDimension>& mesh,
    std::vector<TemplateElement<double, kDimension, kDimension> >&
        template_elements,
    const std::string& solution_filename,
    bool compute_estimator,
    double solver_tolerance = 1.0e-11,
    int maximum_iterations = 500) {
  FEMSpace<double, kDimension> fem_space;
  build_p1_space(mesh, template_elements, fem_space);

  StiffMatrix<kDimension, double> stiffness(fem_space);
  stiffness.algebricAccuracy() = 4;
  stiffness.build();

  FEMFunction<double, kDimension> solution(fem_space);
  Vector<double> right_hand_side;
  Operator::L2Discretize(&heat_source, fem_space, right_hand_side, 5);

  apply_zero_dirichlet(fem_space,
                       stiffness,
                       solution,
                       right_hand_side);

  AMGSolver solver(stiffness);
  solver.solve(solution,
               right_hand_side,
               solver_tolerance,
               maximum_iterations);
  if (!solution_filename.empty()) {
    solution.writeOpenDXData(solution_filename);
  }

  SolveResult result;
  result.degrees_of_freedom = fem_space.n_dof();
  result.sensor_temperature =
      evaluate_sensor_temperature(fem_space, solution);
  result.estimator = 0.0;
  if (compute_estimator) {
    result.eta_squared =
        compute_residual_indicator(mesh,
                                   fem_space,
                                   solution,
                                   result.strong_residual_mean_squared,
                                   result.cell_residual_squared,
                                   result.interior_jump_squared,
                                   result.neumann_boundary_squared);
    result.estimator = std::sqrt(std::accumulate(
        result.eta_squared.begin(), result.eta_squared.end(), 0.0));
  }
  return result;
}

typedef std::pair<long long, long long> CoordinateKey;

CoordinateKey coordinate_key(const Point<kDimension>& point) {
  const double scale = 1.0e12;
  return CoordinateKey(
      static_cast<long long>(std::llround(scale * point[0])),
      static_cast<long long>(std::llround(scale * point[1])));
}

void prolong_p1_function(
    IrregularMesh<kDimension>& coarse_irregular_mesh,
    FEMSpace<double, kDimension>& coarse_space,
    const FEMFunction<double, kDimension>& coarse_function,
    IrregularMesh<kDimension>& fine_irregular_mesh,
    FEMSpace<double, kDimension>& fine_space,
    FEMFunction<double, kDimension>& fine_function) {
  IrregularMeshPair<kDimension> mesh_pair(
      coarse_irregular_mesh, fine_irregular_mesh);
  ActiveElementPairIterator<kDimension> pair =
      mesh_pair.beginActiveElementPair();
  ActiveElementPairIterator<kDimension> end =
      mesh_pair.endActiveElementPair();
  for (; pair != end; ++pair) {
    const HElement<kDimension>& coarse_h_element = pair(0);
    const HElement<kDimension>& fine_h_element = pair(1);
    if (pair.state() != ActiveElementPairIterator<kDimension>::GREAT_THAN &&
        pair.state() != ActiveElementPairIterator<kDimension>::EQUAL) {
      throw std::runtime_error(
          "The DWR enrichment mesh must refine the primal mesh.");
    }
    const Element<double, kDimension>& coarse_element =
        coarse_space.element(coarse_h_element.index);
    const Element<double, kDimension>& fine_element =
        fine_space.element(fine_h_element.index);
    for (unsigned int local_dof = 0;
         local_dof < fine_element.dof().size(); ++local_dof) {
      const int fine_dof = fine_element.dof()[local_dof];
      const Point<kDimension>& point =
          fine_space.dofInfo(fine_dof).interp_point;
      fine_function(fine_dof) =
          coarse_function.value(point, coarse_element);
    }
  }
}

DWRResult compute_dwr_indicator(
    IrregularMesh<kDimension>& coarse_irregular_mesh,
    std::vector<TemplateElement<double, kDimension, kDimension> >&
        template_elements,
    const std::string& primal_solution_filename) {
  RegularMesh<kDimension>& coarse_mesh =
      coarse_irregular_mesh.regularMesh();

  // Coarse primal: R_h(U_h) = A_h U_h - b_h = 0 in the P1 test space.
  FEMSpace<double, kDimension> coarse_space;
  build_p1_space(coarse_mesh, template_elements, coarse_space);
  StiffMatrix<kDimension, double> coarse_stiffness(coarse_space);
  coarse_stiffness.algebricAccuracy() = 4;
  coarse_stiffness.build();
  Vector<double> coarse_right_hand_side;
  Operator::L2Discretize(
      &heat_source, coarse_space, coarse_right_hand_side, 5);
  FEMFunction<double, kDimension> coarse_primal(coarse_space);
  apply_zero_dirichlet(coarse_space,
                       coarse_stiffness,
                       coarse_primal,
                       coarse_right_hand_side);
  AMGSolver coarse_solver(coarse_stiffness);
  coarse_solver.solve(
      coarse_primal, coarse_right_hand_side, 1.0e-11, 500);
  coarse_primal.writeOpenDXData(primal_solution_filename);
  std::vector<double> coarse_strong_residual_mean_squared;
  std::vector<double> coarse_cell_residual_squared;
  std::vector<double> coarse_interior_jump_squared;
  std::vector<double> coarse_neumann_boundary_squared;
  const std::vector<double> coarse_residual_eta_squared =
      compute_residual_indicator(coarse_mesh,
                                 coarse_space,
                                 coarse_primal,
                                 coarse_strong_residual_mean_squared,
                                 coarse_cell_residual_squared,
                                 coarse_interior_jump_squared,
                                 coarse_neumann_boundary_squared);

  // Enrich by one uniform h-refinement and keep P1.  This produces a fully
  // discrete fine-space dual while retaining AFEPack's nested geometry tree.
  IrregularMesh<kDimension> fine_irregular_mesh(coarse_irregular_mesh);
  fine_irregular_mesh.globalRefine(1);
  fine_irregular_mesh.semiregularize();
  fine_irregular_mesh.regularize(false);
  RegularMesh<kDimension>& fine_mesh = fine_irregular_mesh.regularMesh();
  FEMSpace<double, kDimension> fine_space;
  build_p1_space(fine_mesh, template_elements, fine_space);

  StiffMatrix<kDimension, double> fine_jacobian(fine_space);
  fine_jacobian.algebricAccuracy() = 4;
  fine_jacobian.build();
  Vector<double> fine_dual_right_hand_side;
  Operator::L2Discretize(&sensor_functional_gradient,
                         fine_space,
                         fine_dual_right_hand_side,
                         5);
  for (unsigned int dof = 0;
       dof < fine_dual_right_hand_side.size(); ++dof) {
    fine_dual_right_hand_side(dof) *= -1.0;
  }
  FEMFunction<double, kDimension> fine_dual(fine_space);
  apply_zero_dirichlet(fine_space,
                       fine_jacobian,
                       fine_dual,
                       fine_dual_right_hand_side);
  AMGSolver fine_solver(fine_jacobian);
  fine_solver.solve(fine_dual, fine_dual_right_hand_side, 1.0e-11, 500);

  FEMFunction<double, kDimension> fine_primal(fine_space);
  prolong_p1_function(coarse_irregular_mesh,
                      coarse_space,
                      coarse_primal,
                      fine_irregular_mesh,
                      fine_space,
                      fine_primal);

  // Pi_h Z_H is obtained at the nested coarse nodes and then prolonged back
  // to V_H.  The correction Z_H - I_h^H Pi_h Z_H is the part that the primal
  // test space cannot represent.
  std::map<CoordinateKey, double> fine_nodal_dual;
  for (unsigned int dof = 0; dof < fine_space.n_dof(); ++dof) {
    fine_nodal_dual[coordinate_key(fine_space.dofInfo(dof).interp_point)] =
        fine_dual(dof);
  }
  FEMFunction<double, kDimension> coarse_dual_projection(coarse_space);
  for (unsigned int dof = 0; dof < coarse_space.n_dof(); ++dof) {
    const CoordinateKey key =
        coordinate_key(coarse_space.dofInfo(dof).interp_point);
    const std::map<CoordinateKey, double>::const_iterator found =
        fine_nodal_dual.find(key);
    if (found == fine_nodal_dual.end()) {
      throw std::runtime_error(
          "Cannot match a coarse node with the enriched dual mesh.");
    }
    coarse_dual_projection(dof) = found->second;
  }
  FEMFunction<double, kDimension> fine_dual_projection(fine_space);
  prolong_p1_function(coarse_irregular_mesh,
                      coarse_space,
                      coarse_dual_projection,
                      fine_irregular_mesh,
                      fine_space,
                      fine_dual_projection);

  DWRResult result;
  const int coarse_element_count = coarse_mesh.n_geometry(kDimension);
  result.signed_element_contribution.assign(coarse_element_count, 0.0);
  result.absolute_element_indicator.assign(coarse_element_count, 0.0);
  result.primal_residual_norm.assign(coarse_element_count, 0.0);
  result.dual_solution_rms.assign(coarse_element_count, 0.0);
  result.dual_correction_rms.assign(coarse_element_count, 0.0);
  result.residual_dual_correction_product.assign(coarse_element_count, 0.0);
  result.residual_absolute_dual_product.assign(coarse_element_count, 0.0);
  result.signed_estimate = 0.0;
  result.absolute_sum = 0.0;
  result.sensor_temperature =
      evaluate_sensor_temperature(coarse_space, coarse_primal);
  result.primal_degrees_of_freedom = coarse_space.n_dof();
  result.enriched_dual_degrees_of_freedom = fine_space.n_dof();

  IrregularMeshPair<kDimension> mesh_pair(
      coarse_irregular_mesh, fine_irregular_mesh);
  ActiveElementPairIterator<kDimension> pair =
      mesh_pair.beginActiveElementPair();
  ActiveElementPairIterator<kDimension> end =
      mesh_pair.endActiveElementPair();
  const int fine_element_count = fine_mesh.n_geometry(kDimension);
  const int fine_edge_count = fine_mesh.n_geometry(kDimension - 1);
  std::vector<int> fine_to_coarse(fine_element_count, -1);
  for (; pair != end; ++pair) {
    const int coarse_element_index = pair(0).index;
    const int fine_element_index = pair(1).index;
    fine_to_coarse[fine_element_index] = coarse_element_index;
  }

  // Localize the weak residual by integration by parts.  Using raw element
  // weak forms before this step produces large artificial cancellations
  // across faces.  The signed sum remains
  // R_H(I_h^H U_h)(Z_H-I_h^H Pi_h Z_H).  Marking uses the absolute value of
  // each localized contribution.  Residual/dual norm products are exported
  // only as diagnostics.
  std::vector<std::vector<int> > edge_elements(fine_edge_count);
  std::vector<double> dual_correction_l2_squared(
      coarse_element_count, 0.0);
  std::vector<double> dual_solution_l2_squared(
      coarse_element_count, 0.0);
  std::vector<double> accumulated_volume(coarse_element_count, 0.0);
  for (int fine_element_index = 0;
       fine_element_index < fine_element_count; ++fine_element_index) {
    const int coarse_element_index = fine_to_coarse[fine_element_index];
    if (coarse_element_index < 0) {
      throw std::runtime_error(
          "Cannot map an enriched element to its coarse parent.");
    }
    const Element<double, kDimension>& fine_element =
        fine_space.element(fine_element_index);
    const GeometryBM& geometry =
        fine_mesh.geometry(kDimension, fine_element_index);
    const QuadratureInfo<kDimension>& quadrature =
        fine_element.findQuadratureInfo(5);
    const std::vector<Point<kDimension> > points =
        fine_element.local_to_global(quadrature.quadraturePoint());
    const std::vector<double> jacobians =
        fine_element.local_to_global_jacobian(
            quadrature.quadraturePoint());
    const double reference_volume = fine_element.templateElement().volume();
    for (int q = 0; q < quadrature.n_quadraturePoint(); ++q) {
      const double dual_correction =
          fine_dual.value(points[q], fine_element) -
          fine_dual_projection.value(points[q], fine_element);
      const double dual_value = fine_dual.value(points[q], fine_element);
      const double weight = quadrature.weight(q) * jacobians[q] *
                            reference_volume;
      result.signed_element_contribution[coarse_element_index] -=
          weight * heat_source(points[q]) * dual_correction;
      dual_correction_l2_squared[coarse_element_index] +=
          weight * dual_correction * dual_correction;
      dual_solution_l2_squared[coarse_element_index] +=
          weight * dual_value * dual_value;
      accumulated_volume[coarse_element_index] += weight;
    }

    // Twin-triangle macro elements contain an internal diagonal that is not
    // listed in RegularMesh::geometry(1).  Add its flux jump explicitly.
    if (geometry.n_vertex() == 4) {
      const Point<kDimension>& vertex_0 =
          fine_mesh.point(geometry.vertex(0));
      const Point<kDimension>& vertex_1 =
          fine_mesh.point(geometry.vertex(1));
      const Point<kDimension>& vertex_2 =
          fine_mesh.point(geometry.vertex(2));
      const Point<kDimension>& vertex_3 =
          fine_mesh.point(geometry.vertex(3));
      Point<kDimension> upper_centroid;
      Point<kDimension> lower_centroid;
      Point<kDimension> diagonal_midpoint;
      for (int component = 0; component < kDimension; ++component) {
        upper_centroid[component] =
            (vertex_0[component] + vertex_1[component] +
             vertex_2[component]) / 3.0;
        lower_centroid[component] =
            (vertex_0[component] + vertex_2[component] +
             vertex_3[component]) / 3.0;
        diagonal_midpoint[component] =
            0.5 * (vertex_0[component] + vertex_2[component]);
      }
      const double diagonal_length =
          std::sqrt(distance_squared(vertex_0, vertex_2));
      std::vector<double> upper_outward_normal(2);
      upper_outward_normal[0] =
          (vertex_2[1] - vertex_0[1]) / diagonal_length;
      upper_outward_normal[1] =
          -(vertex_2[0] - vertex_0[0]) / diagonal_length;
      const double toward_upper =
          upper_outward_normal[0] *
              (upper_centroid[0] - diagonal_midpoint[0]) +
          upper_outward_normal[1] *
              (upper_centroid[1] - diagonal_midpoint[1]);
      if (toward_upper > 0.0) {
        upper_outward_normal[0] *= -1.0;
        upper_outward_normal[1] *= -1.0;
      }
      const std::vector<double> gradient_upper =
          fine_primal.gradient(upper_centroid, fine_element);
      const std::vector<double> gradient_lower =
          fine_primal.gradient(lower_centroid, fine_element);
      const double flux_jump = kConductivity *
          ((gradient_upper[0] - gradient_lower[0]) *
               upper_outward_normal[0] +
           (gradient_upper[1] - gradient_lower[1]) *
               upper_outward_normal[1]);
      const double dual_correction =
          fine_dual.value(diagonal_midpoint, fine_element) -
          fine_dual_projection.value(diagonal_midpoint, fine_element);
      result.signed_element_contribution[coarse_element_index] +=
          diagonal_length * flux_jump * dual_correction;
    }

    for (int local_edge = 0;
         local_edge < geometry.n_boundary(); ++local_edge) {
      edge_elements[geometry.boundary(local_edge)].push_back(
          fine_element_index);
    }
  }

  for (int edge_index = 0; edge_index < fine_edge_count; ++edge_index) {
    const GeometryBM& edge =
        fine_mesh.geometry(kDimension - 1, edge_index);
    const std::vector<int>& neighbors = edge_elements[edge_index];
    const Point<kDimension> midpoint = edge_midpoint(fine_mesh, edge);
    const double length = edge_length(fine_mesh, edge);
    if (neighbors.size() == 2) {
      const int first_index = neighbors[0];
      const int second_index = neighbors[1];
      const GeometryBM& first_geometry =
          fine_mesh.geometry(kDimension, first_index);
      const GeometryBM& second_geometry =
          fine_mesh.geometry(kDimension, second_index);
      const std::vector<double> first_normal = outward_unit_normal(
          fine_mesh, edge, first_geometry);
      const std::vector<double> second_normal = outward_unit_normal(
          fine_mesh, edge, second_geometry);
      const std::vector<double> first_gradient = fine_primal.gradient(
          midpoint, fine_space.element(first_index));
      const std::vector<double> second_gradient = fine_primal.gradient(
          midpoint, fine_space.element(second_index));
      const double flux_jump = kConductivity *
          (first_gradient[0] * first_normal[0] +
           first_gradient[1] * first_normal[1] +
           second_gradient[0] * second_normal[0] +
           second_gradient[1] * second_normal[1]);
      const double dual_correction =
          fine_dual.value(midpoint, fine_space.element(first_index)) -
          fine_dual_projection.value(
              midpoint, fine_space.element(first_index));
      const double edge_contribution =
          length * flux_jump * dual_correction;
      result.signed_element_contribution[fine_to_coarse[first_index]] +=
          0.5 * edge_contribution;
      result.signed_element_contribution[fine_to_coarse[second_index]] +=
          0.5 * edge_contribution;
    } else if (neighbors.size() == 1) {
      if (edge.boundaryMark() != 0) {
        continue;
      }
      const int element_index = neighbors[0];
      const GeometryBM& geometry =
          fine_mesh.geometry(kDimension, element_index);
      const std::vector<double> normal = outward_unit_normal(
          fine_mesh, edge, geometry);
      const std::vector<double> gradient = fine_primal.gradient(
          midpoint, fine_space.element(element_index));
      const double normal_flux = kConductivity *
          (gradient[0] * normal[0] + gradient[1] * normal[1]);
      const double dual_correction =
          fine_dual.value(midpoint, fine_space.element(element_index)) -
          fine_dual_projection.value(
              midpoint, fine_space.element(element_index));
      result.signed_element_contribution[
          fine_to_coarse[element_index]] +=
          length * normal_flux * dual_correction;
    } else {
      throw std::runtime_error(
          "Unexpected enriched-edge adjacency in the DWR localization.");
    }
  }

  for (int element_index = 0;
       element_index < coarse_element_count; ++element_index) {
    const double contribution =
        result.signed_element_contribution[element_index];
    result.absolute_element_indicator[element_index] =
        std::abs(contribution);
    // These factors and products are diagnostic.  The actual DWR marking
    // score is abs(contribution), which retains the localized volume and
    // face residual action assembled above.
    result.primal_residual_norm[element_index] =
        std::sqrt(coarse_cell_residual_squared[element_index]);
    result.dual_correction_rms[element_index] =
        accumulated_volume[element_index] > 0.0
            ? std::sqrt(dual_correction_l2_squared[element_index] /
                        accumulated_volume[element_index])
            : 0.0;
    result.dual_solution_rms[element_index] =
        accumulated_volume[element_index] > 0.0
            ? std::sqrt(dual_solution_l2_squared[element_index] /
                        accumulated_volume[element_index])
            : 0.0;
    result.residual_dual_correction_product[element_index] =
        result.primal_residual_norm[element_index] *
        result.dual_correction_rms[element_index];
    result.residual_absolute_dual_product[element_index] =
        result.primal_residual_norm[element_index] *
        result.dual_solution_rms[element_index];
    result.signed_estimate += contribution;
    result.absolute_sum += std::abs(contribution);
  }
  return result;
}

void write_dwr_indicator_data(
    RegularMesh<kDimension>& mesh,
    const DWRResult& dwr,
    const std::vector<bool>& marked,
    const std::string& filename) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  output << "# element centroid_x centroid_y signed_dwr abs_dwr "
            "primal_residual_norm dual_correction_rms "
            "dual_solution_rms residual_dual_correction_product "
            "residual_absolute_dual_product marked\n"
         << std::scientific << std::setprecision(12);
  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& geometry =
        mesh.geometry(kDimension, element_index);
    Point<kDimension> centroid;
    centroid[0] = 0.0;
    centroid[1] = 0.0;
    for (int vertex = 0; vertex < geometry.n_vertex(); ++vertex) {
      centroid[0] += mesh.point(geometry.vertex(vertex))[0];
      centroid[1] += mesh.point(geometry.vertex(vertex))[1];
    }
    centroid[0] /= geometry.n_vertex();
    centroid[1] /= geometry.n_vertex();
    output << element_index << ' ' << centroid[0] << ' ' << centroid[1]
           << ' ' << dwr.signed_element_contribution[element_index]
           << ' ' << dwr.absolute_element_indicator[element_index]
           << ' ' << dwr.primal_residual_norm[element_index]
           << ' ' << dwr.dual_correction_rms[element_index]
           << ' ' << dwr.dual_solution_rms[element_index]
           << ' ' << dwr.residual_dual_correction_product[element_index]
           << ' ' << dwr.residual_absolute_dual_product[element_index]
           << ' ' << (marked[element_index] ? 1 : 0) << '\n';
  }
}

void write_functional_comparison_data(
    const std::vector<FunctionalRecord>& residual_records,
    const std::vector<FunctionalRecord>& dual_records,
    const std::vector<FunctionalRecord>& dwr_records,
    double reference_value,
    double reference_change) {
  const std::string filename = "summary/functional_comparison.dat";
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  output << "# spectral_reference " << std::scientific
         << std::setprecision(14) << reference_value << '\n'
         << "# spectral_reference_change " << reference_change << '\n'
         << "# strategy level elements dofs functional absolute_error\n";
  const auto write_records = [&](const char* strategy,
                                 const std::vector<FunctionalRecord>& records) {
    for (const FunctionalRecord& record : records) {
      output << strategy << ' ' << record.level << ' '
             << record.elements << ' ' << record.degrees_of_freedom << ' '
             << record.value << ' '
             << std::abs(record.value - reference_value) << '\n';
    }
  };
  write_records("residual", residual_records);
  write_records("dual", dual_records);
  write_records("dwr", dwr_records);
}


}  // namespace

int main(int argc, char* argv[]) {
  const std::string root_mesh = argc >= 2 ? argv[1] : "D";
  const bool uniform_reference_mode =
      argc == 4 && std::string(argv[2]) == "--uniform-reference";
  const bool custom_rounds_mode =
      argc == 4 && std::string(argv[2]) == "--rounds";
  const bool comparison_rounds_mode =
      argc == 5 && std::string(argv[2]) == "--comparison-rounds";
  const bool three_way_comparison_rounds_mode =
      argc == 6 && std::string(argv[2]) == "--comparison-rounds";

  if (argc != 2 && !uniform_reference_mode && !custom_rounds_mode &&
      !comparison_rounds_mode && !three_way_comparison_rounds_mode) {
    std::cerr
        << "usage: " << argv[0]
        << " ROOT_MESH\n"
        << "       " << argv[0]
        << " ROOT_MESH --rounds N\n"
        << "       " << argv[0]
        << " ROOT_MESH --comparison-rounds RESIDUAL_ROUNDS DWR_ROUNDS\n"
        << "       " << argv[0]
        << " ROOT_MESH --comparison-rounds RESIDUAL_ROUNDS "
           "DUAL_ROUNDS DWR_ROUNDS\n"
        << "       " << argv[0]
        << " ROOT_MESH --uniform-reference LEVEL\n";
    return EXIT_FAILURE;
  }
  int uniform_reference_level = 0;
  int adaptation_rounds = kDefaultAdaptationRounds;
  int residual_rounds = adaptation_rounds;
  int dual_rounds = adaptation_rounds;
  int dwr_rounds = adaptation_rounds;
  if (uniform_reference_mode) {
    uniform_reference_level = std::stoi(argv[3]);
    if (uniform_reference_level < 0) {
      std::cerr << "uniform reference level must be nonnegative\n";
      return EXIT_FAILURE;
    }
  }
  if (custom_rounds_mode) {
    adaptation_rounds = std::stoi(argv[3]);
    if (adaptation_rounds <= 0) {
      std::cerr << "the number of adaptation rounds must be positive\n";
      return EXIT_FAILURE;
    }
    residual_rounds = adaptation_rounds;
    dual_rounds = adaptation_rounds;
    dwr_rounds = adaptation_rounds;
  }
  if (comparison_rounds_mode) {
    residual_rounds = std::stoi(argv[3]);
    dwr_rounds = std::stoi(argv[4]);
    dual_rounds = dwr_rounds;
    if (residual_rounds <= 0 || dwr_rounds <= 0) {
      std::cerr << "the comparison round counts must be positive\n";
      return EXIT_FAILURE;
    }
  }
  if (three_way_comparison_rounds_mode) {
    residual_rounds = std::stoi(argv[3]);
    dual_rounds = std::stoi(argv[4]);
    dwr_rounds = std::stoi(argv[5]);
    if (residual_rounds <= 0 || dual_rounds <= 0 || dwr_rounds <= 0) {
      std::cerr << "the comparison round counts must be positive\n";
      return EXIT_FAILURE;
    }
  }

  try {
    prepare_output_directories();

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

    if (uniform_reference_mode) {
      HGeometryTree<kDimension> reference_hierarchy;
      reference_hierarchy.readEasyMesh(root_mesh);
      IrregularMesh<kDimension> reference_irregular_mesh(
          reference_hierarchy);
      reference_irregular_mesh.globalRefine(uniform_reference_level);
      reference_irregular_mesh.semiregularize();
      reference_irregular_mesh.regularize(false);
      RegularMesh<kDimension>& reference_mesh =
          reference_irregular_mesh.regularMesh();
      const SolveResult reference_solution = solve_temperature(
          reference_mesh,
          template_elements,
          "",
          false,
          1.0e-9,
          500);
      const double spectral_reference =
          spectral_reference_sensor_temperature(
              kReferenceModes, kReferenceIntervals);
      const std::string output_filename =
          "summary/uniform_reference_level_" +
          std::to_string(uniform_reference_level) + ".dat";
      std::ofstream reference_output(output_filename.c_str());
      reference_output << std::scientific << std::setprecision(14)
                       << "# uniform refinement reference\n"
                       << "level " << uniform_reference_level << '\n'
                       << "elements "
                       << reference_mesh.n_geometry(kDimension) << '\n'
                       << "dofs "
                       << reference_solution.degrees_of_freedom << '\n'
                       << "functional "
                       << reference_solution.sensor_temperature << '\n'
                       << "spectral_reference "
                       << spectral_reference << '\n'
                       << "absolute_difference "
                       << std::abs(reference_solution.sensor_temperature -
                                   spectral_reference)
                       << '\n';
      std::cout << std::setprecision(14)
                << "Uniform reference solve\n"
                << "  refinement level    : "
                << uniform_reference_level << '\n'
                << "  elements            : "
                << reference_mesh.n_geometry(kDimension) << '\n'
                << "  dofs                : "
                << reference_solution.degrees_of_freedom << '\n'
                << "  J_h                 : "
                << reference_solution.sensor_temperature << '\n'
                << "  spectral J_ref      : "
                << spectral_reference << '\n'
                << "  absolute difference : "
                << std::abs(reference_solution.sensor_temperature -
                            spectral_reference)
                << '\n'
                << "  output              : "
                << output_filename << '\n';
      return EXIT_SUCCESS;
    }

    HGeometryTree<kDimension> hierarchy;
    hierarchy.readEasyMesh(root_mesh);
    IrregularMesh<kDimension> irregular_mesh(hierarchy);
    irregular_mesh.semiregularize();
    irregular_mesh.regularize(false);

    std::vector<AdaptationReport> reports;
    for (int round = 1; round <= residual_rounds; ++round) {
      RegularMesh<kDimension>& current_mesh = irregular_mesh.regularMesh();
      const std::string state_suffix = std::to_string(round - 1);
      const std::string mesh_prefix =
          "meshes/residual/level_" + state_suffix;
      const std::string temperature_filename =
          "fields/residual/temperature_level_" + state_suffix + ".dx";

      current_mesh.writeOpenDXData(mesh_prefix + ".dx");
      if (round == 1) {
        const std::vector<double> source_samples =
            sample_at_element_centroids(current_mesh,
                                        &heat_source_at_point);
        const std::vector<double> sensor_samples =
            sample_at_element_centroids(current_mesh, &sensor_weight);
        write_cell_field_data(current_mesh,
                              source_samples,
                              "fields/problem/source_profile.dat");
        write_cell_field_data(current_mesh,
                              sensor_samples,
                              "fields/problem/sensor_weight.dat");
      }

      const SolveResult current = solve_temperature(
          current_mesh,
          template_elements,
          temperature_filename,
          true);
      const int elements_before =
          current_mesh.n_geometry(kDimension);
      // For this introductory example, refine from the strong element
      // residual only.  The full estimator is still computed and exported,
      // but its flux-jump and boundary terms are deliberately diagnostic:
      // they are face terms and can form rings or isolated marked cells.
      const std::vector<double>& marking_indicator =
          current.cell_residual_squared;
      const ProportionalMarking marking =
          mark_largest_fraction(marking_indicator, kRefineFraction);
      const std::vector<bool>& marked = marking.marked;
      const int marked_count = static_cast<int>(
          std::count(marked.begin(), marked.end(), true));
      const std::string indicator_stem =
          "fields/residual/indicator_round_" + std::to_string(round);

      write_indicator_data(current_mesh,
                           current.eta_squared,
                           current.strong_residual_mean_squared,
                           current.cell_residual_squared,
                           current.interior_jump_squared,
                           current.neumann_boundary_squared,
                           marked,
                           indicator_stem + ".dat");
      if (round == residual_rounds) {
        write_indicator_data(current_mesh,
                             current.eta_squared,
                             current.strong_residual_mean_squared,
                             current.cell_residual_squared,
                             current.interior_jump_squared,
                             current.neumann_boundary_squared,
                             marked,
                             "fields/residual/indicator_final.dat");
      }

      adapt_with_afepack_quantile(irregular_mesh,
                                  current_mesh,
                                  marking_indicator,
                                  marking.tolerance);

      irregular_mesh.semiregularize();
      irregular_mesh.regularize(false);
      RegularMesh<kDimension>& adapted_mesh = irregular_mesh.regularMesh();
      reports.push_back(AdaptationReport{
          round,
          elements_before,
          current.degrees_of_freedom,
          current.estimator,
          current.sensor_temperature,
          marked_count,
          static_cast<int>(adapted_mesh.n_geometry(kDimension)),
          count_twin_triangles(adapted_mesh)});
    }

    RegularMesh<kDimension>& mesh_after = irregular_mesh.regularMesh();
    mesh_after.writeOpenDXData("meshes/residual/final.dx");
    const SolveResult after = solve_temperature(
        mesh_after,
        template_elements,
        "fields/residual/temperature_final.dx",
        false);

    const double reference_value = spectral_reference_sensor_temperature(
        kReferenceModes, kReferenceIntervals);
    const double reference_check = spectral_reference_sensor_temperature(
        120, kReferenceIntervals);
    const double reference_change =
        std::abs(reference_value - reference_check);
    std::vector<FunctionalRecord> residual_records;
    for (const AdaptationReport& report : reports) {
      residual_records.push_back(FunctionalRecord{
          report.round - 1,
          report.elements_before,
          report.dofs_before,
          report.sensor_temperature});
    }
    residual_records.push_back(FunctionalRecord{
        residual_rounds,
        static_cast<int>(mesh_after.n_geometry(kDimension)),
        after.degrees_of_freedom,
        after.sensor_temperature});

    const std::string residual_history_filename =
        "summary/residual_history.dat";
    std::ofstream functional_history(residual_history_filename.c_str());
    if (!functional_history) {
      throw std::runtime_error("Cannot open " + residual_history_filename);
    }
    functional_history
        << "# spectral_reference " << reference_value << '\n'
        << "# refinement_level elements dofs sensor_temperature absolute_error\n"
        << std::scientific << std::setprecision(12);
    for (const FunctionalRecord& record : residual_records) {
      functional_history << record.level << ' '
                         << record.elements << ' '
                         << record.degrees_of_freedom << ' '
                         << record.value << ' '
                         << std::abs(record.value - reference_value) << '\n';
    }

    // The dual uses its own hierarchy, initialized again from the root mesh.
    // It is therefore never solved on the residual-adapted primal mesh.
    HGeometryTree<kDimension> dual_hierarchy;
    dual_hierarchy.readEasyMesh(root_mesh);
    IrregularMesh<kDimension> dual_irregular_mesh(dual_hierarchy);
    dual_irregular_mesh.semiregularize();
    dual_irregular_mesh.regularize(false);

    std::vector<DualAdaptationReport> dual_reports;
    for (int round = 1; round <= dual_rounds; ++round) {
      RegularMesh<kDimension>& dual_mesh =
          dual_irregular_mesh.regularMesh();
      const std::string round_string = std::to_string(round);
      dual_mesh.writeOpenDXData(
          "meshes/dual/level_" + std::to_string(round - 1) + ".dx");

      const SolveResult primal_on_dual_mesh = solve_temperature(
          dual_mesh,
          template_elements,
          "fields/dual/temperature_round_" + round_string + ".dx",
          false);
      const DiscreteDualResult dual_result = solve_discrete_dual(
          dual_mesh,
          template_elements,
          "fields/dual/signed_round_" + round_string);
      const ProportionalMarking dual_marking = mark_largest_fraction(
          dual_result.element_rms_magnitude, kRefineFraction);
      const std::vector<bool>& dual_marked = dual_marking.marked;
      const int dual_marked_count = static_cast<int>(
          std::count(dual_marked.begin(), dual_marked.end(), true));
      write_dual_magnitude_data(
          dual_mesh,
          dual_result.element_rms_magnitude,
          dual_marked,
          "fields/dual/magnitude_round_" + round_string + ".dat");

      adapt_with_afepack_quantile(dual_irregular_mesh,
                                  dual_mesh,
                                  dual_result.element_rms_magnitude,
                                  dual_marking.tolerance);

      dual_irregular_mesh.semiregularize();
      dual_irregular_mesh.regularize(false);
      dual_reports.push_back(DualAdaptationReport{
          round,
          dual_result.elements,
          dual_result.degrees_of_freedom,
          primal_on_dual_mesh.sensor_temperature,
          dual_result.minimum,
          dual_result.maximum,
          dual_marked_count,
          static_cast<int>(
              dual_irregular_mesh.regularMesh().n_geometry(kDimension))});
    }

    RegularMesh<kDimension>& dual_mesh_after =
        dual_irregular_mesh.regularMesh();
    dual_mesh_after.writeOpenDXData("meshes/dual/final.dx");
    const DiscreteDualResult dual_after = solve_discrete_dual(
        dual_mesh_after,
        template_elements,
        "fields/dual/signed_final");
    const SolveResult primal_on_dual_mesh_after = solve_temperature(
        dual_mesh_after,
        template_elements,
        "fields/dual/temperature_final.dx",
        false);
    const std::vector<bool> no_dual_marks(
        dual_after.element_rms_magnitude.size(), false);
    write_dual_magnitude_data(dual_mesh_after,
                              dual_after.element_rms_magnitude,
                              no_dual_marks,
                              "fields/dual/magnitude_final.dat");

    std::vector<FunctionalRecord> dual_records;
    for (const DualAdaptationReport& report : dual_reports) {
      dual_records.push_back(FunctionalRecord{
          report.round - 1,
          report.elements_before,
          report.dofs_before,
          report.sensor_temperature});
    }
    dual_records.push_back(FunctionalRecord{
        dual_rounds,
        static_cast<int>(dual_mesh_after.n_geometry(kDimension)),
        primal_on_dual_mesh_after.degrees_of_freedom,
        primal_on_dual_mesh_after.sensor_temperature});

    const std::string dual_history_filename = "summary/dual_history.dat";
    std::ofstream dual_history(dual_history_filename.c_str());
    if (!dual_history) {
      throw std::runtime_error("Cannot open " + dual_history_filename);
    }
    dual_history
        << "# level elements dofs functional dual_min dual_max "
           "marked elements_after\n"
        << std::scientific << std::setprecision(12);
    for (const DualAdaptationReport& report : dual_reports) {
      dual_history << report.round - 1 << ' '
                   << report.elements_before << ' '
                   << report.dofs_before << ' '
                   << report.sensor_temperature << ' '
                   << report.minimum << ' '
                   << report.maximum << ' '
                   << report.marked_elements << ' '
                   << report.elements_after << '\n';
    }
    dual_history << dual_rounds << ' '
                 << dual_mesh_after.n_geometry(kDimension) << ' '
                 << primal_on_dual_mesh_after.degrees_of_freedom << ' '
                 << primal_on_dual_mesh_after.sensor_temperature
                 << " nan nan 0 "
                 << dual_mesh_after.n_geometry(kDimension) << '\n';

    // A third independent sequence uses an absolute adjoint-weighted residual
    // for marking.  The primal remains P1 on the current mesh, while the
    // fully discrete dual is solved in P1 on one uniformly h-refined copy.
    // This keeps the adapted primal DOF counts comparable with the two
    // proportional-marking sequences above.
    HGeometryTree<kDimension> dwr_hierarchy;
    dwr_hierarchy.readEasyMesh(root_mesh);
    IrregularMesh<kDimension> dwr_irregular_mesh(dwr_hierarchy);
    dwr_irregular_mesh.semiregularize();
    dwr_irregular_mesh.regularize(false);

    std::vector<DWRAdaptationReport> dwr_reports;
    for (int round = 1; round <= dwr_rounds; ++round) {
      RegularMesh<kDimension>& dwr_mesh =
          dwr_irregular_mesh.regularMesh();
      const int dwr_elements_before =
          dwr_mesh.n_geometry(kDimension);
      const std::string round_string = std::to_string(round);
      dwr_mesh.writeOpenDXData(
          "meshes/dwr/level_" + std::to_string(round - 1) + ".dx");
      const DWRResult dwr = compute_dwr_indicator(
          dwr_irregular_mesh,
          template_elements,
          "fields/dwr/temperature_round_" + round_string + ".dx");
      const std::vector<double>& dwr_marking_indicator =
          dwr.absolute_element_indicator;
      const ProportionalMarking dwr_marking = mark_largest_fraction(
          dwr_marking_indicator, kRefineFraction);
      const std::vector<bool>& dwr_marked = dwr_marking.marked;
      const int dwr_marked_count = static_cast<int>(
          std::count(dwr_marked.begin(), dwr_marked.end(), true));
      write_dwr_indicator_data(
          dwr_mesh,
          dwr,
          dwr_marked,
          "fields/dwr/indicator_round_" + round_string + ".dat");
      if (round == dwr_rounds) {
        write_dwr_indicator_data(
            dwr_mesh,
            dwr,
            dwr_marked,
            "fields/dwr/indicator_final.dat");
      }

      adapt_with_afepack_quantile(dwr_irregular_mesh,
                                  dwr_mesh,
                                  dwr_marking_indicator,
                                  dwr_marking.tolerance);

      dwr_irregular_mesh.semiregularize();
      dwr_irregular_mesh.regularize(false);
      dwr_reports.push_back(DWRAdaptationReport{
          round,
          dwr_elements_before,
          dwr.primal_degrees_of_freedom,
          dwr.sensor_temperature,
          dwr.signed_estimate,
          dwr.absolute_sum,
          dwr.enriched_dual_degrees_of_freedom,
          dwr_marked_count,
          static_cast<int>(
              dwr_irregular_mesh.regularMesh().n_geometry(kDimension))});
    }

    RegularMesh<kDimension>& dwr_mesh_after =
        dwr_irregular_mesh.regularMesh();
    dwr_mesh_after.writeOpenDXData("meshes/dwr/final.dx");
    const SolveResult primal_on_dwr_mesh_after = solve_temperature(
        dwr_mesh_after,
        template_elements,
        "fields/dwr/temperature_final.dx",
        false);
    std::vector<FunctionalRecord> dwr_records;
    for (const DWRAdaptationReport& report : dwr_reports) {
      dwr_records.push_back(FunctionalRecord{
          report.round - 1,
          report.elements_before,
          report.dofs_before,
          report.sensor_temperature});
    }
    dwr_records.push_back(FunctionalRecord{
        dwr_rounds,
        static_cast<int>(dwr_mesh_after.n_geometry(kDimension)),
        primal_on_dwr_mesh_after.degrees_of_freedom,
        primal_on_dwr_mesh_after.sensor_temperature});

    const std::string dwr_history_filename = "summary/dwr_history.dat";
    std::ofstream dwr_history(dwr_history_filename.c_str());
    if (!dwr_history) {
      throw std::runtime_error("Cannot open " + dwr_history_filename);
    }
    dwr_history
        << "# spectral_reference " << std::scientific
        << std::setprecision(14) << reference_value << '\n'
        << "# level elements dofs functional signed_estimate "
           "corrected_functional absolute_error corrected_absolute_error\n";
    for (const DWRAdaptationReport& report : dwr_reports) {
      const double corrected_functional =
          report.sensor_temperature + report.signed_estimate;
      dwr_history << report.round - 1 << ' '
                  << report.elements_before << ' '
                  << report.dofs_before << ' '
                  << report.sensor_temperature << ' '
                  << report.signed_estimate << ' '
                  << corrected_functional << ' '
                  << std::abs(report.sensor_temperature - reference_value)
                  << ' '
                  << std::abs(corrected_functional - reference_value)
                  << '\n';
    }
    dwr_history << dwr_rounds << ' '
                << dwr_mesh_after.n_geometry(kDimension) << ' '
                << primal_on_dwr_mesh_after.degrees_of_freedom << ' '
                << primal_on_dwr_mesh_after.sensor_temperature
                << " nan nan "
                << std::abs(primal_on_dwr_mesh_after.sensor_temperature -
                            reference_value)
                << " nan\n";
    write_functional_comparison_data(residual_records,
                                     dual_records,
                                     dwr_records,
                                     reference_value,
                                     reference_change);

    std::cout << "\nRun summary\n"
              << "  residual / dual / DWR rounds : "
              << residual_rounds << " / " << dual_rounds << " / "
              << dwr_rounds << '\n'
              << "  marked fraction              : "
              << 100.0 * kRefineFraction << "%\n"
              << "  spectral reference J_ref     : "
              << std::scientific << std::setprecision(8)
              << reference_value << "\n\n"
              << std::left << std::setw(12) << "strategy"
              << std::right << std::setw(7) << "level"
              << std::setw(11) << "elements"
              << std::setw(9) << "dofs"
              << std::setw(16) << "J(T_h)"
              << std::setw(16) << "|J-J_ref|" << '\n';
    const auto print_records = [&](const char* strategy,
                                   const std::vector<FunctionalRecord>& records) {
      for (const FunctionalRecord& record : records) {
        std::cout << std::left << std::setw(12) << strategy
                  << std::right << std::setw(7) << record.level
                  << std::setw(11) << record.elements
                  << std::setw(9) << record.degrees_of_freedom
                  << std::setw(16) << record.value
                  << std::setw(16)
                  << std::abs(record.value - reference_value) << '\n';
      }
    };
    print_records("residual", residual_records);
    print_records("dual", dual_records);
    print_records("DWR", dwr_records);

    std::cout << "\nDWR correction diagnostic\n"
              << std::left << std::setw(8) << "level"
              << std::right << std::setw(16) << "eta_signed"
              << std::setw(16) << "J+eta"
              << std::setw(16) << "|J+eta-Jref|" << '\n';
    for (const DWRAdaptationReport& report : dwr_reports) {
      const double corrected_functional =
          report.sensor_temperature + report.signed_estimate;
      std::cout << std::left << std::setw(8) << report.round - 1
                << std::right << std::setw(16) << report.signed_estimate
                << std::setw(16) << corrected_functional
                << std::setw(16)
                << std::abs(corrected_functional - reference_value)
                << '\n';
    }
  } catch (const std::exception& error) {
    std::cerr << "two_heater_goal_adaptive: " << error.what() << '\n';
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
