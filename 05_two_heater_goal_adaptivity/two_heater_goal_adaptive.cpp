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
constexpr double kDualColorGamma = 4.0;
constexpr int kReferenceModes = 160;
constexpr int kReferenceIntervals = 4096;

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

struct Bounds {
  double xmin;
  double xmax;
  double ymin;
  double ymax;
};

Bounds mesh_bounds(const RegularMesh<kDimension>& mesh) {
  Bounds bounds{mesh.point(0)[0], mesh.point(0)[0],
                mesh.point(0)[1], mesh.point(0)[1]};
  for (unsigned int i = 1; i < mesh.n_point(); ++i) {
    bounds.xmin = std::min(bounds.xmin, mesh.point(i)[0]);
    bounds.xmax = std::max(bounds.xmax, mesh.point(i)[0]);
    bounds.ymin = std::min(bounds.ymin, mesh.point(i)[1]);
    bounds.ymax = std::max(bounds.ymax, mesh.point(i)[1]);
  }
  return bounds;
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

double svg_x(double x, const Bounds& bounds) {
  return 55.0 + 690.0 * (x - bounds.xmin) /
                    (bounds.xmax - bounds.xmin);
}

double svg_y(double y, const Bounds& bounds) {
  return 745.0 - 690.0 * (y - bounds.ymin) /
                     (bounds.ymax - bounds.ymin);
}

void write_mesh_edges_svg(const RegularMesh<kDimension>& mesh,
                          const std::string& filename,
                          const std::string& caption) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  const Bounds bounds = mesh_bounds(mesh);
  output << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 800 800\">\n"
         << "<rect width=\"800\" height=\"800\" fill=\"white\"/>\n"
         << "<text x=\"400\" y=\"30\" text-anchor=\"middle\" "
         << "font-family=\"sans-serif\" font-size=\"20\">"
         << caption << "</text>\n";
  output << std::setprecision(10);
  for (unsigned int edge_index = 0;
       edge_index < mesh.n_geometry(kDimension - 1); ++edge_index) {
    const GeometryBM& edge =
        mesh.geometry(kDimension - 1, edge_index);
    const Point<kDimension>& first = mesh.point(edge.vertex(0));
    const Point<kDimension>& second = mesh.point(edge.vertex(1));
    output << "<line x1=\"" << svg_x(first[0], bounds)
           << "\" y1=\"" << svg_y(first[1], bounds)
           << "\" x2=\"" << svg_x(second[0], bounds)
           << "\" y2=\"" << svg_y(second[1], bounds)
           << "\" stroke=\"#40566b\" stroke-width=\"0.75\"/>\n";
  }
  output << "</svg>\n";
}

void write_clean_mesh_svg(const RegularMesh<kDimension>& mesh,
                          const std::string& filename) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  const Bounds bounds = mesh_bounds(mesh);
  output << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 800 800\">\n"
         << "<rect width=\"800\" height=\"800\" fill=\"white\"/>\n";
  output << std::setprecision(10);
  for (unsigned int edge_index = 0;
       edge_index < mesh.n_geometry(kDimension - 1); ++edge_index) {
    const GeometryBM& edge = mesh.geometry(kDimension - 1, edge_index);
    const Point<kDimension>& first = mesh.point(edge.vertex(0));
    const Point<kDimension>& second = mesh.point(edge.vertex(1));
    output << "<line x1=\"" << svg_x(first[0], bounds)
           << "\" y1=\"" << svg_y(first[1], bounds)
           << "\" x2=\"" << svg_x(second[0], bounds)
           << "\" y2=\"" << svg_y(second[1], bounds)
           << "\" stroke=\"#52677a\" stroke-width=\"0.58\"/>\n";
  }
  output << "</svg>\n";
}

void write_indicator_svg(const RegularMesh<kDimension>& mesh,
                         const std::vector<double>& eta_squared,
                         const std::vector<bool>& marked,
                         const std::string& filename,
                         const std::string& caption,
                         bool show_marked = true) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  const Bounds bounds = mesh_bounds(mesh);
  const double maximum =
      *std::max_element(eta_squared.begin(), eta_squared.end());

  output << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 800 800\">\n"
         << "<rect width=\"800\" height=\"800\" fill=\"white\"/>\n"
         << "<text x=\"400\" y=\"30\" text-anchor=\"middle\" "
         << "font-family=\"sans-serif\" font-size=\"20\">"
         << caption << "</text>\n";
  output << std::setprecision(10);

  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& element = mesh.geometry(kDimension, element_index);
    const double ratio = maximum > 0.0
                             ? std::sqrt(eta_squared[element_index] / maximum)
                             : 0.0;
    const int red = 245;
    const int green = static_cast<int>(242.0 - 150.0 * ratio);
    const int blue = static_cast<int>(224.0 - 185.0 * ratio);
    output << "<polygon points=\"";
    for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
      const Point<kDimension>& point = mesh.point(element.vertex(vertex));
      output << svg_x(point[0], bounds) << ',' << svg_y(point[1], bounds)
             << ' ';
    }
    const bool is_marked = show_marked && marked[element_index];
    output << "\" fill=\"rgb(" << red << ',' << green << ',' << blue
           << ")\" stroke=\""
           << (is_marked ? "#9d1825" : "#738496")
           << "\" stroke-width=\""
           << (is_marked ? "1.8" : "0.55") << "\"/>\n";
  }
  if (show_marked) {
    output << "<rect x=\"570\" y=\"760\" width=\"18\" height=\"12\" "
           << "fill=\"#f55c45\" stroke=\"#9d1825\"/>\n"
           << "<text x=\"596\" y=\"771\" font-family=\"sans-serif\" "
           << "font-size=\"14\">marked elements</text>\n";
  }
  output << "</svg>\n";
}

enum class CleanIndicatorPalette {
  residual,
  dual,
  dwr
};

void write_clean_indicator_svg(
    const RegularMesh<kDimension>& mesh,
    const std::vector<double>& indicator,
    const std::vector<bool>& marked,
    const std::string& filename,
    CleanIndicatorPalette palette,
    double display_gamma) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  const Bounds bounds = mesh_bounds(mesh);
  const double maximum =
      *std::max_element(indicator.begin(), indicator.end());

  int dark_red = 217;
  int dark_green = 93;
  int dark_blue = 32;
  if (palette == CleanIndicatorPalette::dual) {
    dark_red = 15;
    dark_green = 118;
    dark_blue = 110;
  } else if (palette == CleanIndicatorPalette::dwr) {
    dark_red = 109;
    dark_green = 66;
    dark_blue = 194;
  }
  const int light_red = 249;
  const int light_green = 248;
  const int light_blue = 241;

  output << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 800 800\">\n"
         << "<rect width=\"800\" height=\"800\" fill=\"white\"/>\n";
  output << std::setprecision(10);
  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& element = mesh.geometry(kDimension, element_index);
    const double normalized = maximum > 0.0
                                  ? indicator[element_index] / maximum
                                  : 0.0;
    const double focused = std::pow(normalized, display_gamma);
    const int red = static_cast<int>(std::lround(
        light_red + (dark_red - light_red) * focused));
    const int green = static_cast<int>(std::lround(
        light_green + (dark_green - light_green) * focused));
    const int blue = static_cast<int>(std::lround(
        light_blue + (dark_blue - light_blue) * focused));
    output << "<polygon points=\"";
    for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
      const Point<kDimension>& point = mesh.point(element.vertex(vertex));
      output << svg_x(point[0], bounds) << ',' << svg_y(point[1], bounds)
             << ' ';
    }
    const bool is_marked = marked[element_index];
    output << "\" fill=\"rgb(" << red << ',' << green << ',' << blue
           << ")\" stroke=\""
           << (is_marked ? "#172033" : "#8796a5")
           << "\" stroke-width=\""
           << (is_marked ? "2.0" : "0.48")
           << "\" stroke-linejoin=\"round\"/>\n";
  }
  output << "</svg>\n";
}

void write_indicator_data(const RegularMesh<kDimension>& mesh,
                          const std::vector<double>& eta_squared,
                          const std::vector<double>& cell_residual_squared,
                          const std::vector<double>& interior_jump_squared,
                          const std::vector<double>& neumann_boundary_squared,
                          const std::vector<bool>& marked,
                          const std::string& filename) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  output << "# element centroid_x centroid_y eta_squared cell_squared "
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
           << cell_residual_squared[element_index] << ' '
           << interior_jump_squared[element_index] << ' '
           << neumann_boundary_squared[element_index] << ' '
           << (marked[element_index] ? 1 : 0) << '\n';
  }
}

void write_dual_magnitude_svg(
    const RegularMesh<kDimension>& mesh,
    const std::vector<double>& magnitude,
    const std::vector<bool>& marked,
    const std::string& filename,
    const std::string& caption,
    bool show_marked) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }
  const Bounds bounds = mesh_bounds(mesh);
  const double maximum =
      *std::max_element(magnitude.begin(), magnitude.end());

  output << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 800 800\">\n"
         << "<rect width=\"800\" height=\"800\" fill=\"white\"/>\n"
         << "<text x=\"400\" y=\"30\" text-anchor=\"middle\" "
         << "font-family=\"sans-serif\" font-size=\"20\">"
         << caption << "</text>\n";
  output << std::setprecision(10);

  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& element = mesh.geometry(kDimension, element_index);
    const double normalized =
        maximum > 0.0 ? magnitude[element_index] / maximum : 0.0;
    // The display exponent suppresses the small global tail of the elliptic
    // dual while leaving the computed values unchanged.
    const double focused = std::pow(normalized, kDualColorGamma);
    const int red = static_cast<int>(244.0 - 210.0 * focused);
    const int green = static_cast<int>(247.0 - 130.0 * focused);
    const int blue = static_cast<int>(242.0 - 100.0 * focused);
    output << "<polygon points=\"";
    for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
      const Point<kDimension>& point = mesh.point(element.vertex(vertex));
      output << svg_x(point[0], bounds) << ',' << svg_y(point[1], bounds)
             << ' ';
    }
    const bool is_marked = show_marked && marked[element_index];
    output << "\" fill=\"rgb(" << red << ',' << green << ',' << blue
           << ")\" stroke=\""
           << (is_marked ? "#9d1825" : "#738496")
           << "\" stroke-width=\""
           << (is_marked ? "1.8" : "0.45") << "\"/>\n";
  }

  const double distractor_x = svg_x(kDistractorCenterX, bounds);
  const double distractor_y = svg_y(kDistractorCenterY, bounds);
  const double distractor_radius_x =
      690.0 * kDistractorSigmaX / (bounds.xmax - bounds.xmin);
  const double distractor_radius_y =
      690.0 * kDistractorSigmaY / (bounds.ymax - bounds.ymin);
  const double target_x = svg_x(kTargetCenterX, bounds);
  const double target_y = svg_y(kTargetCenterY, bounds);
  const double target_radius =
      690.0 * kTargetSigma / (bounds.xmax - bounds.xmin);
  const double sensor_x = svg_x(kSensorCenterX, bounds);
  const double sensor_y = svg_y(kSensorCenterY, bounds);
  const double sensor_radius =
      690.0 * kSensorSigma / (bounds.xmax - bounds.xmin);
  output << "<ellipse cx=\"" << distractor_x << "\" cy=\"" << distractor_y
         << "\" rx=\"" << distractor_radius_x
         << "\" ry=\"" << distractor_radius_y
         << "\" fill=\"none\" stroke=\"#d97706\" stroke-width=\"2.5\"/>\n"
         << "<text x=\"" << distractor_x + distractor_radius_x + 7.0
         << "\" y=\"" << distractor_y - 5.0
         << "\" font-family=\"sans-serif\" font-size=\"13\" "
         << "fill=\"#92400e\">strong distractor heater</text>\n"
         << "<circle cx=\"" << target_x << "\" cy=\"" << target_y
         << "\" r=\"" << target_radius
         << "\" fill=\"none\" stroke=\"#ca8a04\" stroke-width=\"2\" "
         << "stroke-dasharray=\"6 4\"/>\n"
         << "<text x=\"" << target_x - target_radius - 8.0
         << "\" y=\"" << target_y + target_radius + 16.0
         << "\" text-anchor=\"end\" font-family=\"sans-serif\" "
         << "font-size=\"13\" fill=\"#854d0e\">smooth target heater</text>\n"
         << "<circle cx=\"" << sensor_x << "\" cy=\"" << sensor_y
         << "\" r=\"" << sensor_radius
         << "\" fill=\"none\" stroke=\"#116466\" stroke-width=\"3\"/>\n"
         << "<text x=\"" << sensor_x - sensor_radius - 8.0
         << "\" y=\"" << sensor_y - 5.0
         << "\" text-anchor=\"end\" font-family=\"sans-serif\" "
         << "font-size=\"13\" fill=\"#116466\">sensor</text>\n";
  if (show_marked) {
    output << "<rect x=\"505\" y=\"760\" width=\"18\" height=\"12\" "
           << "fill=\"#f55c45\" stroke=\"#9d1825\"/>\n"
           << "<text x=\"531\" y=\"771\" font-family=\"sans-serif\" "
           << "font-size=\"13\">marked elements</text>\n";
  } else {
    output << "<text x=\"55\" y=\"775\" font-family=\"sans-serif\" "
           << "font-size=\"13\">color contrast exponent = "
           << kDualColorGamma << "; max |psi| = " << maximum
           << "</text>\n";
  }
  output << "</svg>\n";
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

void write_discrete_dual_svg(
    const RegularMesh<kDimension>& mesh,
    FEMSpace<double, kDimension>& fem_space,
    const FEMFunction<double, kDimension>& dual,
    const std::string& filename) {
  std::ofstream output(filename.c_str());
  if (!output) {
    throw std::runtime_error("Cannot open " + filename);
  }

  std::vector<double> values(mesh.n_geometry(kDimension), 0.0);
  double minimum = 0.0;
  double maximum = 0.0;
  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& geometry =
        mesh.geometry(kDimension, element_index);
    Point<kDimension> centroid;
    centroid[0] = 0.0;
    centroid[1] = 0.0;
    for (int vertex = 0; vertex < geometry.n_vertex(); ++vertex) {
      const Point<kDimension>& point = mesh.point(geometry.vertex(vertex));
      centroid[0] += point[0];
      centroid[1] += point[1];
    }
    centroid[0] /= geometry.n_vertex();
    centroid[1] /= geometry.n_vertex();
    values[element_index] =
        dual.value(centroid, fem_space.element(element_index));
    minimum = std::min(minimum, values[element_index]);
    maximum = std::max(maximum, values[element_index]);
  }

  const double scale = std::max(std::abs(minimum), std::abs(maximum));
  const Bounds bounds = mesh_bounds(mesh);
  output << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 800 800\">\n"
         << "<rect width=\"800\" height=\"800\" fill=\"white\"/>\n"
         << "<text x=\"400\" y=\"30\" text-anchor=\"middle\" "
         << "font-family=\"sans-serif\" font-size=\"20\">"
         << "Fully discrete dual: A^T psi_h = -grad J_h</text>\n";
  output << std::setprecision(10);

  for (unsigned int element_index = 0;
       element_index < mesh.n_geometry(kDimension); ++element_index) {
    const GeometryBM& element = mesh.geometry(kDimension, element_index);
    const double signed_ratio =
        scale > 0.0 ? values[element_index] / scale : 0.0;
    const double magnitude = std::abs(signed_ratio);
    int red = 248;
    int green = 247;
    int blue = 242;
    if (signed_ratio < 0.0) {
      red = static_cast<int>(248.0 - 60.0 * magnitude);
      green = static_cast<int>(242.0 - 175.0 * magnitude);
      blue = static_cast<int>(232.0 - 190.0 * magnitude);
    } else {
      red = static_cast<int>(235.0 - 175.0 * magnitude);
      green = static_cast<int>(245.0 - 105.0 * magnitude);
      blue = static_cast<int>(250.0 - 25.0 * magnitude);
    }
    output << "<polygon points=\"";
    for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
      const Point<kDimension>& point = mesh.point(element.vertex(vertex));
      output << svg_x(point[0], bounds) << ',' << svg_y(point[1], bounds)
             << ' ';
    }
    output << "\" fill=\"rgb(" << red << ',' << green << ',' << blue
           << ")\" stroke=\"#738496\" stroke-width=\"0.45\"/>\n";
  }

  const double sensor_x = svg_x(kSensorCenterX, bounds);
  const double sensor_y = svg_y(kSensorCenterY, bounds);
  const double sensor_radius =
      690.0 * kSensorSigma / (bounds.xmax - bounds.xmin);
  output << "<circle cx=\"" << sensor_x << "\" cy=\"" << sensor_y
         << "\" r=\"" << sensor_radius
         << "\" fill=\"none\" stroke=\"#116466\" stroke-width=\"3\"/>\n"
         << "<text x=\"" << sensor_x + sensor_radius + 8.0
         << "\" y=\"" << sensor_y - 5.0
         << "\" font-family=\"sans-serif\" font-size=\"14\" "
         << "fill=\"#116466\">sensor footprint</text>\n"
         << "<text x=\"55\" y=\"775\" font-family=\"sans-serif\" "
         << "font-size=\"14\">min psi_h = " << minimum
         << ", max psi_h = " << maximum << "</text>\n"
         << "</svg>\n";
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
  write_discrete_dual_svg(mesh,
                          fem_space,
                          dual,
                          output_stem + ".svg");
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
  std::ofstream output("functional_comparison.dat");
  if (!output) {
    throw std::runtime_error("Cannot open functional_comparison.dat");
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

double plot_x(int dofs, int minimum_dofs, int maximum_dofs) {
  if (maximum_dofs == minimum_dofs) {
    return 80.0;
  }
  return 80.0 + 650.0 * (dofs - minimum_dofs) /
                    static_cast<double>(maximum_dofs - minimum_dofs);
}

void write_functional_value_svg(
    const std::vector<FunctionalRecord>& residual_records,
    const std::vector<FunctionalRecord>& dual_records,
    const std::vector<FunctionalRecord>& dwr_records,
    double reference_value) {
  int minimum_dofs = residual_records.front().degrees_of_freedom;
  int maximum_dofs = minimum_dofs;
  double minimum_value = residual_records.front().value;
  double maximum_value = minimum_value;
  for (const std::vector<FunctionalRecord>* records :
       {&residual_records, &dual_records, &dwr_records}) {
    for (const FunctionalRecord& record : *records) {
      minimum_dofs = std::min(minimum_dofs, record.degrees_of_freedom);
      maximum_dofs = std::max(maximum_dofs, record.degrees_of_freedom);
      minimum_value = std::min(minimum_value, record.value);
      maximum_value = std::max(maximum_value, record.value);
    }
  }
  const double value_span = std::max(maximum_value - minimum_value, 1.0e-12);
  const double padding = 0.15 * value_span;
  minimum_value -= padding;
  maximum_value += padding;
  const auto y_coordinate = [&](double value) {
    return 440.0 - 350.0 * (value - minimum_value) /
                       (maximum_value - minimum_value);
  };

  std::ofstream output("functional_value_vs_dofs.svg");
  if (!output) {
    throw std::runtime_error("Cannot open functional_value_vs_dofs.svg");
  }
  output << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 800 520\">\n"
         << "<rect width=\"800\" height=\"520\" fill=\"white\"/>\n"
         << "<text x=\"400\" y=\"32\" text-anchor=\"middle\" "
         << "font-family=\"sans-serif\" font-size=\"21\">"
         << "Target functional versus degrees of freedom</text>\n"
         << "<line x1=\"80\" y1=\"440\" x2=\"730\" y2=\"440\" "
         << "stroke=\"#374151\"/>\n"
         << "<line x1=\"80\" y1=\"90\" x2=\"80\" y2=\"440\" "
         << "stroke=\"#374151\"/>\n";
  output << std::setprecision(10);
  for (int tick = 0; tick <= 5; ++tick) {
    const double ratio = tick / 5.0;
    const int dofs = static_cast<int>(std::lround(
        minimum_dofs + ratio * (maximum_dofs - minimum_dofs)));
    const double x = plot_x(dofs, minimum_dofs, maximum_dofs);
    const double value = minimum_value + ratio *
                         (maximum_value - minimum_value);
    const double y = y_coordinate(value);
    output << "<line x1=\"" << x << "\" y1=\"440\" x2=\"" << x
           << "\" y2=\"446\" stroke=\"#374151\"/>\n"
           << "<text x=\"" << x << "\" y=\"466\" text-anchor=\"middle\" "
           << "font-family=\"sans-serif\" font-size=\"12\">" << dofs
           << "</text>\n"
           << "<line x1=\"74\" y1=\"" << y << "\" x2=\"80\" y2=\""
           << y << "\" stroke=\"#374151\"/>\n"
           << "<text x=\"68\" y=\"" << y + 4.0
           << "\" text-anchor=\"end\" font-family=\"sans-serif\" "
           << "font-size=\"12\">" << std::fixed << std::setprecision(5)
           << value << "</text>\n";
  }
  output << "<text x=\"725\" y=\"84\" text-anchor=\"end\" "
         << "font-family=\"sans-serif\" font-size=\"11\">"
         << "spectral reference J = " << std::fixed << std::setprecision(8)
         << reference_value << " (outside vertical scale)</text>\n";

  const auto write_series = [&](const std::vector<FunctionalRecord>& records,
                                const char* color,
                                int marker) {
    for (std::size_t index = 1; index < records.size(); ++index) {
      const FunctionalRecord& previous = records[index - 1];
      const FunctionalRecord& current = records[index];
      output << "<line x1=\""
             << plot_x(previous.degrees_of_freedom,
                       minimum_dofs,
                       maximum_dofs)
             << "\" y1=\"" << y_coordinate(previous.value)
             << "\" x2=\""
             << plot_x(current.degrees_of_freedom,
                       minimum_dofs,
                       maximum_dofs)
             << "\" y2=\"" << y_coordinate(current.value)
             << "\" stroke=\"" << color
             << "\" stroke-width=\"3\"/>\n";
    }
    for (const FunctionalRecord& record : records) {
      const double x = plot_x(record.degrees_of_freedom,
                              minimum_dofs,
                              maximum_dofs);
      const double y = y_coordinate(record.value);
      if (marker == 0) {
        output << "<circle cx=\"" << x << "\" cy=\"" << y
               << "\" r=\"5\" fill=\"" << color << "\"/>\n";
      } else if (marker == 1) {
        output << "<rect x=\"" << x - 4.5 << "\" y=\"" << y - 4.5
               << "\" width=\"9\" height=\"9\" fill=\"" << color
               << "\"/>\n";
      } else {
        output << "<polygon points=\"" << x << ',' << y - 5.5 << ' '
               << x - 5.0 << ',' << y + 4.5 << ' '
               << x + 5.0 << ',' << y + 4.5
               << "\" fill=\"" << color << "\"/>\n";
      }
    }
  };
  write_series(residual_records, "#d97706", 0);
  write_series(dual_records, "#0f766e", 1);
  write_series(dwr_records, "#7c3aed", 2);
  output << "<circle cx=\"115\" cy=\"62\" r=\"5\" fill=\"#d97706\"/>\n"
         << "<text x=\"128\" y=\"66\" font-family=\"sans-serif\" "
         << "font-size=\"13\">residual refinement</text>\n"
         << "<rect x=\"295\" y=\"57\" width=\"10\" height=\"10\" "
         << "fill=\"#0f766e\"/>\n"
         << "<text x=\"313\" y=\"66\" font-family=\"sans-serif\" "
         << "font-size=\"13\">dual-magnitude refinement</text>\n"
         << "<polygon points=\"585,56 580,66 590,66\" "
         << "fill=\"#7c3aed\"/>\n"
         << "<text x=\"598\" y=\"66\" font-family=\"sans-serif\" "
         << "font-size=\"13\">DWR refinement</text>\n"
         << "<text x=\"405\" y=\"500\" text-anchor=\"middle\" "
         << "font-family=\"sans-serif\" font-size=\"14\">degrees of freedom</text>\n"
         << "<text x=\"18\" y=\"270\" text-anchor=\"middle\" "
         << "font-family=\"sans-serif\" font-size=\"14\" "
         << "transform=\"rotate(-90 18 270)\">J(T_h)</text>\n"
         << "</svg>\n";
}

void write_functional_error_svg(
    const std::vector<FunctionalRecord>& residual_records,
    const std::vector<FunctionalRecord>& dual_records,
    const std::vector<FunctionalRecord>& dwr_records,
    double reference_value) {
  int minimum_dofs = residual_records.front().degrees_of_freedom;
  int maximum_dofs = minimum_dofs;
  double minimum_error =
      std::abs(residual_records.front().value - reference_value);
  double maximum_error = minimum_error;
  for (const std::vector<FunctionalRecord>* records :
       {&residual_records, &dual_records, &dwr_records}) {
    for (const FunctionalRecord& record : *records) {
      minimum_dofs = std::min(minimum_dofs, record.degrees_of_freedom);
      maximum_dofs = std::max(maximum_dofs, record.degrees_of_freedom);
      const double error = std::abs(record.value - reference_value);
      minimum_error = std::min(minimum_error, error);
      maximum_error = std::max(maximum_error, error);
    }
  }
  double minimum_log_error = std::log10(minimum_error);
  double maximum_log_error = std::log10(maximum_error);
  const double log_span =
      std::max(maximum_log_error - minimum_log_error, 1.0);
  minimum_log_error -= 0.08 * log_span;
  maximum_log_error += 0.08 * log_span;
  const auto y_coordinate = [&](double error) {
    return 440.0 - 350.0 * (std::log10(error) - minimum_log_error) /
                       (maximum_log_error - minimum_log_error);
  };

  std::ofstream output("functional_error_vs_dofs.svg");
  if (!output) {
    throw std::runtime_error("Cannot open functional_error_vs_dofs.svg");
  }
  output << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 800 520\">\n"
         << "<rect width=\"800\" height=\"520\" fill=\"white\"/>\n"
         << "<text x=\"400\" y=\"32\" text-anchor=\"middle\" "
         << "font-family=\"sans-serif\" font-size=\"21\">"
         << "Target functional error versus degrees of freedom (log scale)"
         << "</text>\n"
         << "<line x1=\"80\" y1=\"440\" x2=\"730\" y2=\"440\" "
         << "stroke=\"#374151\"/>\n"
         << "<line x1=\"80\" y1=\"90\" x2=\"80\" y2=\"440\" "
         << "stroke=\"#374151\"/>\n";
  output << std::setprecision(10);
  for (int tick = 0; tick <= 5; ++tick) {
    const double ratio = tick / 5.0;
    const int dofs = static_cast<int>(std::lround(
        minimum_dofs + ratio * (maximum_dofs - minimum_dofs)));
    const double x = plot_x(dofs, minimum_dofs, maximum_dofs);
    output << "<line x1=\"" << x << "\" y1=\"440\" x2=\"" << x
           << "\" y2=\"446\" stroke=\"#374151\"/>\n"
           << "<text x=\"" << x << "\" y=\"466\" text-anchor=\"middle\" "
           << "font-family=\"sans-serif\" font-size=\"12\">" << dofs
           << "</text>\n";
  }
  for (int tick = 0; tick <= 5; ++tick) {
    const double ratio = tick / 5.0;
    const double error = std::pow(
        10.0,
        minimum_log_error + ratio *
            (maximum_log_error - minimum_log_error));
    const double y = y_coordinate(error);
    output << "<line x1=\"80\" y1=\"" << y << "\" x2=\"730\" y2=\""
           << y << "\" stroke=\"#d1d5db\" stroke-width=\"0.8\"/>\n"
           << "<text x=\"68\" y=\"" << y + 4.0
           << "\" text-anchor=\"end\" font-family=\"sans-serif\" "
           << "font-size=\"12\">" << std::scientific
           << std::setprecision(1)
           << error << "</text>\n";
  }
  output << std::defaultfloat << std::setprecision(10);
  const auto write_series = [&](const std::vector<FunctionalRecord>& records,
                                const char* color,
                                int marker) {
    for (std::size_t index = 1; index < records.size(); ++index) {
      const FunctionalRecord& previous = records[index - 1];
      const FunctionalRecord& current = records[index];
      output << "<line x1=\""
             << plot_x(previous.degrees_of_freedom,
                       minimum_dofs,
                       maximum_dofs)
             << "\" y1=\""
             << y_coordinate(std::abs(previous.value - reference_value))
             << "\" x2=\""
             << plot_x(current.degrees_of_freedom,
                       minimum_dofs,
                       maximum_dofs)
             << "\" y2=\""
             << y_coordinate(std::abs(current.value - reference_value))
             << "\" stroke=\"" << color
             << "\" stroke-width=\"3\"/>\n";
    }
    for (const FunctionalRecord& record : records) {
      const double x = plot_x(record.degrees_of_freedom,
                              minimum_dofs,
                              maximum_dofs);
      const double y = y_coordinate(std::abs(record.value - reference_value));
      if (marker == 0) {
        output << "<circle cx=\"" << x << "\" cy=\"" << y
               << "\" r=\"5\" fill=\"" << color << "\"/>\n";
      } else if (marker == 1) {
        output << "<rect x=\"" << x - 4.5 << "\" y=\"" << y - 4.5
               << "\" width=\"9\" height=\"9\" fill=\"" << color
               << "\"/>\n";
      } else {
        output << "<polygon points=\"" << x << ',' << y - 5.5 << ' '
               << x - 5.0 << ',' << y + 4.5 << ' '
               << x + 5.0 << ',' << y + 4.5
               << "\" fill=\"" << color << "\"/>\n";
      }
    }
  };
  write_series(residual_records, "#d97706", 0);
  write_series(dual_records, "#0f766e", 1);
  write_series(dwr_records, "#7c3aed", 2);
  output << "<circle cx=\"115\" cy=\"62\" r=\"5\" fill=\"#d97706\"/>\n"
         << "<text x=\"128\" y=\"66\" font-family=\"sans-serif\" "
         << "font-size=\"13\">residual refinement</text>\n"
         << "<rect x=\"295\" y=\"57\" width=\"10\" height=\"10\" "
         << "fill=\"#0f766e\"/>\n"
         << "<text x=\"313\" y=\"66\" font-family=\"sans-serif\" "
         << "font-size=\"13\">dual-magnitude refinement</text>\n"
         << "<polygon points=\"585,56 580,66 590,66\" "
         << "fill=\"#7c3aed\"/>\n"
         << "<text x=\"598\" y=\"66\" font-family=\"sans-serif\" "
         << "font-size=\"13\">DWR refinement</text>\n"
         << "<text x=\"405\" y=\"500\" text-anchor=\"middle\" "
         << "font-family=\"sans-serif\" font-size=\"14\">degrees of freedom</text>\n"
         << "<text x=\"18\" y=\"270\" text-anchor=\"middle\" "
         << "font-family=\"sans-serif\" font-size=\"14\" "
         << "transform=\"rotate(-90 18 270)\">|J-J_ref|</text>\n"
         << "</svg>\n";
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
          "uniform_reference_level_" +
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

    std::filesystem::create_directories("teaching");

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
          round == 1 ? "mesh_before" : "mesh_round_" + state_suffix;
      const std::string temperature_filename =
          round == 1 ? "temperature_before.dx"
                     : "temperature_round_" + state_suffix + ".dx";

      current_mesh.writeOpenDXData(mesh_prefix + ".dx");
      write_mesh_edges_svg(
          current_mesh,
          mesh_prefix + ".svg",
          "Mesh before residual refinement round " +
              std::to_string(round));
      if (round == 1) {
        const std::vector<bool> no_marks(
            current_mesh.n_geometry(kDimension), false);
        const std::vector<double> source_samples =
            sample_at_element_centroids(current_mesh,
                                        &heat_source_at_point);
        const std::vector<double> sensor_samples =
            sample_at_element_centroids(current_mesh, &sensor_weight);
        write_indicator_svg(
            current_mesh,
            source_samples,
            no_marks,
            "source_profile.svg",
            "PDE source f: strong distractor and smooth target heater",
            false);
        write_indicator_svg(
            current_mesh,
            sensor_samples,
            no_marks,
            "sensor_weight.svg",
            "Normalized sensor weight w",
            false);
        write_clean_indicator_svg(current_mesh,
                                  source_samples,
                                  no_marks,
                                  "teaching/source_profile_clean.svg",
                                  CleanIndicatorPalette::residual,
                                  0.5);
        write_clean_indicator_svg(current_mesh,
                                  sensor_samples,
                                  no_marks,
                                  "teaching/sensor_weight_clean.svg",
                                  CleanIndicatorPalette::dual,
                                  kDualColorGamma);
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
          "residual_indicator_round_" + std::to_string(round);

      write_indicator_data(current_mesh,
                           current.eta_squared,
                           current.cell_residual_squared,
                           current.interior_jump_squared,
                           current.neumann_boundary_squared,
                           marked,
                           indicator_stem + ".dat");
      write_indicator_svg(
          current_mesh,
          marking_indicator,
          marked,
          indicator_stem + ".svg",
          "Cell residual with top 5% marks, round " +
              std::to_string(round));
      if (round == 1 || round == residual_rounds) {
        write_clean_indicator_svg(
            current_mesh,
            marking_indicator,
            marked,
            "teaching/residual_indicator_round_" +
                std::to_string(round) + "_clean.svg",
            CleanIndicatorPalette::residual,
            0.5);
      }
      if (round == residual_rounds) {
        write_indicator_data(current_mesh,
                             current.eta_squared,
                             current.cell_residual_squared,
                             current.interior_jump_squared,
                             current.neumann_boundary_squared,
                             marked,
                             "residual_indicator.dat");
        write_indicator_svg(current_mesh,
                            marking_indicator,
                            marked,
                            "residual_indicator.svg",
                            "Marking indicator: h^2 times cell residual");
        write_indicator_svg(current_mesh,
                            current.strong_residual_mean_squared,
                            marked,
                            "strong_residual.svg",
                            "Strong residual RMS (independent of mesh size)",
                            false);
        write_indicator_svg(current_mesh,
                            current.eta_squared,
                            marked,
                            "residual_raw.svg",
                            "Full residual estimator (diagnostic only)",
                            false);
        write_indicator_svg(current_mesh,
                            current.cell_residual_squared,
                            marked,
                            "residual_cell_component.svg",
                            "Cell residual contribution",
                            false);
        write_indicator_svg(current_mesh,
                            current.interior_jump_squared,
                            marked,
                            "residual_jump_component.svg",
                            "Interior flux-jump contribution",
                            false);
        write_indicator_svg(current_mesh,
                            current.neumann_boundary_squared,
                            marked,
                            "residual_neumann_component.svg",
                            "Neumann-boundary contribution",
                            false);
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
    mesh_after.writeOpenDXData("mesh_after.dx");
    write_mesh_edges_svg(
        mesh_after,
        "mesh_after.svg",
        "Mesh after 5% residual refinement");
    write_clean_mesh_svg(mesh_after, "teaching/residual_mesh_final_clean.svg");
    const SolveResult after = solve_temperature(
        mesh_after,
        template_elements,
        "temperature_after.dx",
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

    std::ofstream functional_history("functional_history.dat");
    if (!functional_history) {
      throw std::runtime_error("Cannot open functional_history.dat");
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
      write_mesh_edges_svg(
          dual_mesh,
          "dual_mesh_round_" + std::to_string(round - 1) + ".svg",
          "Independent dual mesh before magnitude refinement round " +
              round_string);

      const SolveResult primal_on_dual_mesh = solve_temperature(
          dual_mesh,
          template_elements,
          "temperature_on_dual_mesh_round_" + round_string + ".dx",
          false);
      const DiscreteDualResult dual_result = solve_discrete_dual(
          dual_mesh,
          template_elements,
          "dual_signed_round_" + round_string);
      const ProportionalMarking dual_marking = mark_largest_fraction(
          dual_result.element_rms_magnitude, kRefineFraction);
      const std::vector<bool>& dual_marked = dual_marking.marked;
      const int dual_marked_count = static_cast<int>(
          std::count(dual_marked.begin(), dual_marked.end(), true));
      write_dual_magnitude_svg(
          dual_mesh,
          dual_result.element_rms_magnitude,
          dual_marked,
          "dual_magnitude_round_" + round_string + ".svg",
          "Dual magnitude with top 5% marks, round " + round_string,
          true);
      write_dual_magnitude_data(
          dual_mesh,
          dual_result.element_rms_magnitude,
          dual_marked,
          "dual_magnitude_round_" + round_string + ".dat");
      if (round == 1 || round == dual_rounds) {
        write_clean_indicator_svg(
            dual_mesh,
            dual_result.element_rms_magnitude,
            dual_marked,
            "teaching/dual_magnitude_round_" + round_string +
                "_clean.svg",
            CleanIndicatorPalette::dual,
            kDualColorGamma);
      }

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
    dual_mesh_after.writeOpenDXData("dual_mesh_after.dx");
    write_mesh_edges_svg(dual_mesh_after,
                         "dual_mesh_after.svg",
                         "Dual mesh after magnitude-based refinement");
    write_clean_mesh_svg(dual_mesh_after,
                         "teaching/dual_mesh_final_clean.svg");
    const DiscreteDualResult dual_after = solve_discrete_dual(
        dual_mesh_after,
        template_elements,
        "dual_signed_after");
    const SolveResult primal_on_dual_mesh_after = solve_temperature(
        dual_mesh_after,
        template_elements,
        "temperature_on_dual_mesh_after.dx",
        false);
    const std::vector<bool> no_dual_marks(
        dual_after.element_rms_magnitude.size(), false);
    write_dual_magnitude_svg(
        dual_mesh_after,
        dual_after.element_rms_magnitude,
        no_dual_marks,
        "dual_magnitude.svg",
        "Discrete dual magnitude near the diagonal sensor",
        false);
    write_dual_magnitude_data(dual_mesh_after,
                              dual_after.element_rms_magnitude,
                              no_dual_marks,
                              "dual_magnitude.dat");

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
      write_mesh_edges_svg(
          dwr_mesh,
          "dwr_mesh_round_" + std::to_string(round - 1) + ".svg",
          "Mesh before DWR refinement round " + round_string);
      const DWRResult dwr = compute_dwr_indicator(
          dwr_irregular_mesh,
          template_elements,
          "temperature_on_dwr_mesh_round_" + round_string + ".dx");
      const std::vector<double>& dwr_marking_indicator =
          dwr.absolute_element_indicator;
      const ProportionalMarking dwr_marking = mark_largest_fraction(
          dwr_marking_indicator, kRefineFraction);
      const std::vector<bool>& dwr_marked = dwr_marking.marked;
      const int dwr_marked_count = static_cast<int>(
          std::count(dwr_marked.begin(), dwr_marked.end(), true));
      write_indicator_svg(
          dwr_mesh,
          dwr_marking_indicator,
          dwr_marked,
          "dwr_indicator_round_" + round_string + ".svg",
          "Localized DWR contribution, round " +
              round_string);
      write_dwr_indicator_data(
          dwr_mesh,
          dwr,
          dwr_marked,
          "dwr_indicator_round_" + round_string + ".dat");
      if (round == 1 || round == dwr_rounds) {
        write_clean_indicator_svg(
            dwr_mesh,
            dwr_marking_indicator,
            dwr_marked,
            "teaching/dwr_indicator_round_" + round_string +
                "_clean.svg",
            CleanIndicatorPalette::dwr,
            0.5);
      }
      if (round == dwr_rounds) {
        write_indicator_svg(
            dwr_mesh,
            dwr_marking_indicator,
            dwr_marked,
            "dwr_indicator.svg",
            "DWR marking: absolute localized residual action");
        write_indicator_svg(
            dwr_mesh,
            dwr.primal_residual_norm,
            dwr_marked,
            "dwr_primal_residual_factor.svg",
            "DWR factor 1: primal cell residual norm");
        write_indicator_svg(
            dwr_mesh,
            dwr.dual_solution_rms,
            dwr_marked,
            "dwr_dual_weight_factor.svg",
            "DWR factor 2: absolute dual RMS");
        write_indicator_svg(
            dwr_mesh,
            dwr.dual_solution_rms,
            dwr_marked,
            "dwr_dual_solution_diagnostic.svg",
            "Discrete dual RMS (diagnostic)");
        write_dwr_indicator_data(
            dwr_mesh, dwr, dwr_marked, "dwr_indicator.dat");
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
    dwr_mesh_after.writeOpenDXData("dwr_mesh_after.dx");
    write_mesh_edges_svg(dwr_mesh_after,
                         "dwr_mesh_after.svg",
                         "Mesh after DWR refinement");
    write_clean_mesh_svg(dwr_mesh_after,
                         "teaching/dwr_mesh_final_clean.svg");
    const SolveResult primal_on_dwr_mesh_after = solve_temperature(
        dwr_mesh_after,
        template_elements,
        "temperature_on_dwr_mesh_after.dx",
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

    std::ofstream dwr_history("dwr_history.dat");
    if (!dwr_history) {
      throw std::runtime_error("Cannot open dwr_history.dat");
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
    write_functional_value_svg(residual_records,
                               dual_records,
                               dwr_records,
                               reference_value);
    write_functional_error_svg(residual_records,
                               dual_records,
                               dwr_records,
                               reference_value);

    std::cout << std::setprecision(8)
              << "Steady heat-conduction Poisson problem\n"
              << "  distractor source   : (" << kDistractorCenterX << ", "
              << kDistractorCenterY << "), strength "
              << kDistractorStrength << ", sigma ("
              << kDistractorSigmaX << ", "
              << kDistractorSigmaY << ")\n"
              << "  target source       : (" << kTargetCenterX << ", "
              << kTargetCenterY << "), strength "
              << kTargetStrength << ", sigma "
              << kTargetSigma << '\n'
              << "  sensor center       : (" << kSensorCenterX << ", "
              << kSensorCenterY << ")\n"
              << "  sensor sigma        : " << kSensorSigma << '\n'
              << "  marking indicator   : strong cell residual\n"
              << "  refinement fraction : "
              << 100.0 * kRefineFraction
              << "% of active elements per round\n"
              << "  spectral J_ref      : " << reference_value << '\n'
              << "  reference change    : " << reference_change << '\n';
    for (const AdaptationReport& report : reports) {
      std::cout << "  round " << report.round << '\n'
                << "    elements before   : " << report.elements_before
                << '\n'
                << "    dofs before       : " << report.dofs_before << '\n'
                << "    estimator         : " << report.estimator << '\n'
                << "    sensor J(T_h)     : "
                << report.sensor_temperature << '\n'
                << "    |J-J_ref|         : "
                << std::abs(report.sensor_temperature - reference_value)
                << '\n'
                << "    marked elements   : " << report.marked_elements
                << '\n'
                << "    elements after    : " << report.elements_after
                << '\n'
                << "    twin triangles    : "
                << report.twin_triangles_after << '\n';
    }
    std::cout << "  final sensor J(T_h) : " << after.sensor_temperature
              << '\n'
              << "  final |J-J_ref|     : "
              << std::abs(after.sensor_temperature - reference_value)
              << '\n'
              << "  final dofs          : " << after.degrees_of_freedom
              << '\n';
    std::cout << "Independent discrete dual\n"
              << "  marking rule        : same 5% fraction as residual\n"
              << "  color contrast gamma: " << kDualColorGamma << '\n';
    for (const DualAdaptationReport& report : dual_reports) {
      std::cout << "  round " << report.round << '\n'
                << "    elements before   : " << report.elements_before
                << '\n'
                << "    dofs before       : " << report.dofs_before << '\n'
                << "    sensor J(T_h)     : "
                << report.sensor_temperature << '\n'
                << "    |J-J_ref|         : "
                << std::abs(report.sensor_temperature - reference_value)
                << '\n'
                << "    dual range        : [" << report.minimum << ", "
                << report.maximum << "]\n"
                << "    marked elements   : " << report.marked_elements
                << '\n'
                << "    elements after    : " << report.elements_after
                << '\n';
    }
    std::cout << "  final dual range    : [" << dual_after.minimum << ", "
              << dual_after.maximum << "]\n"
              << "  final dual dofs     : "
              << dual_after.degrees_of_freedom << '\n'
              << "  final sensor J(T_h) : "
              << primal_on_dual_mesh_after.sensor_temperature << '\n'
              << "  final |J-J_ref|     : "
              << std::abs(primal_on_dual_mesh_after.sensor_temperature -
                          reference_value)
              << '\n';
    std::cout << "DWR refinement\n"
              << "  enriched dual       : P1 on one uniform h-refinement\n"
              << "  marking rule        : absolute localized DWR contribution\n"
              << "  refinement fraction : "
              << 100.0 * kRefineFraction
              << "% of active elements per round\n";
    for (const DWRAdaptationReport& report : dwr_reports) {
      const double corrected_functional =
          report.sensor_temperature + report.signed_estimate;
      std::cout << "  round " << report.round << '\n'
                << "    elements before   : " << report.elements_before
                << '\n'
                << "    primal dofs       : " << report.dofs_before << '\n'
                << "    enriched dofs     : "
                << report.enriched_dual_dofs << '\n'
                << "    sensor J(T_h)     : "
                << report.sensor_temperature << '\n'
                << "    |J-J_ref|         : "
                << std::abs(report.sensor_temperature - reference_value)
                << '\n'
                << "    signed estimate   : " << report.signed_estimate
                << '\n'
                << "    corrected J       : " << corrected_functional
                << '\n'
                << "    |J_corr-J_ref|    : "
                << std::abs(corrected_functional - reference_value)
                << '\n'
                << "    sum |eta_K|       : " << report.absolute_sum
                << '\n'
                << "    marked elements   : " << report.marked_elements
                << '\n'
                << "    elements after    : " << report.elements_after
                << '\n';
    }
    std::cout << "  final sensor J(T_h) : "
              << primal_on_dwr_mesh_after.sensor_temperature << '\n'
              << "  final |J-J_ref|     : "
              << std::abs(primal_on_dwr_mesh_after.sensor_temperature -
                          reference_value)
              << '\n'
              << "  final primal dofs   : "
              << primal_on_dwr_mesh_after.degrees_of_freedom << '\n';
  } catch (const std::exception& error) {
    std::cerr << "two_heater_goal_adaptive: " << error.what() << '\n';
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
