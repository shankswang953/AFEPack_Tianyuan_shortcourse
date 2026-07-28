#include <AFEPack/AMGSolver.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Operator.h>
#include <AFEPack/TemplateElement.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr int kDimension = 2;
constexpr int kOuterBoundaryMark = 1;
constexpr int kHoleBoundaryMark = 2;
constexpr int kOuterBoundaryPointCount = 128;
constexpr int kHoleBoundaryPointCount = 96;
constexpr double kOuterRadius = 2.0;
constexpr double kHoleRadius = 0.6;
constexpr double kCoefficientLimit = 0.28;
constexpr double kMinimumRadialFactor = 0.42;
constexpr double kMeshSize = 0.10;
constexpr double kPi = 3.14159265358979323846;
constexpr double kTargetHoleArea =
    kPi * kHoleRadius * kHoleRadius;

struct Point2 {
  double x = 0.0;
  double y = 0.0;
};

using Triangle2 = std::array<Point2, 3>;

struct ShapeParameters {
  // r(theta) = scale * (1 + a2 cos(2 theta)
  //                         + b3 sin(3 theta)
  //                         + a4 cos(4 theta)).
  double a2 = 0.20;
  double b3 = -0.15;
  double a4 = 0.10;

  double& operator[](std::size_t index)
  {
    if (index == 0) {
      return a2;
    }
    if (index == 1) {
      return b3;
    }
    return a4;
  }

  double operator[](std::size_t index) const
  {
    if (index == 0) {
      return a2;
    }
    if (index == 1) {
      return b3;
    }
    return a4;
  }
};

struct SolveResult {
  double area = 0.0;
  double temperature_integral = 0.0;
  double mean_temperature = -std::numeric_limits<double>::infinity();
  unsigned int degrees_of_freedom = 0;
  unsigned int elements = 0;
  std::vector<Triangle2> mesh_triangles;
};

struct HistoryRecord {
  int iteration = 0;
  ShapeParameters parameters;
  SolveResult result;
};

struct P1Templates {
  TemplateGeometry<kDimension> geometry;
  CoordTransform<kDimension, kDimension> transform;
  std::unique_ptr<TemplateDOF<kDimension> > dof;
  std::unique_ptr<
      BasisFunctionAdmin<double, kDimension, kDimension> > basis;
  std::vector<
      TemplateElement<double, kDimension, kDimension> > elements;

  P1Templates()
      : elements(1)
  {
    geometry.readData("triangle.tmp_geo");
    transform.readData("triangle.crd_trs");
    dof = std::make_unique<TemplateDOF<kDimension> >(geometry);
    dof->readData("triangle.1.tmp_dof");
    basis = std::make_unique<
        BasisFunctionAdmin<double, kDimension, kDimension> >(*dof);
    basis->readData("triangle.1.bas_fun");
    elements[0].reinit(
        geometry,
        *dof,
        transform,
        *basis);
  }
};

double zero_temperature(const double*)
{
  return 0.0;
}

double unit_heat_source(const double*)
{
  return 1.0;
}

std::string shell_quote(const std::string& text)
{
  std::string quoted = "'";
  for (char character : text) {
    if (character == '\'') {
      quoted += "'\\''";
    } else {
      quoted += character;
    }
  }
  quoted += "'";
  return quoted;
}

fs::path environment_path(
    const char* variable,
    const fs::path& fallback)
{
  const char* value = std::getenv(variable);
  return value != nullptr && *value != '\0' ? fs::path(value) : fallback;
}

double exact_concentric_annulus_mean_temperature()
{
  const double outer_squared = kOuterRadius * kOuterRadius;
  const double hole_squared = kHoleRadius * kHoleRadius;
  const double logarithmic_coefficient =
      (outer_squared - hole_squared)
      / (4.0 * std::log(kOuterRadius / kHoleRadius));
  const double constant =
      0.25 * hole_squared
      - logarithmic_coefficient * std::log(kHoleRadius);
  const auto radial_antiderivative =
      [&](double radius) {
        const double squared = radius * radius;
        return -squared * squared / 16.0
            + logarithmic_coefficient
                * (0.5 * squared * std::log(radius)
                   - 0.25 * squared)
            + 0.5 * constant * squared;
      };
  return 2.0
      * (radial_antiderivative(kOuterRadius)
         - radial_antiderivative(kHoleRadius))
      / (outer_squared - hole_squared);
}

double area_scale(const ShapeParameters& parameters)
{
  const double mode_energy =
      0.5 * (
          parameters.a2 * parameters.a2
          + parameters.b3 * parameters.b3
          + parameters.a4 * parameters.a4);
  // The continuous polar area of the inner hole is fixed.
  return std::sqrt(kTargetHoleArea / kPi)
      / std::sqrt(1.0 + mode_energy);
}

double radial_factor(
    const ShapeParameters& parameters,
    double theta)
{
  return 1.0
      + parameters.a2 * std::cos(2.0 * theta)
      + parameters.b3 * std::sin(3.0 * theta)
      + parameters.a4 * std::cos(4.0 * theta);
}

bool valid_shape(const ShapeParameters& parameters)
{
  for (std::size_t index = 0; index < 3; ++index) {
    if (std::abs(parameters[index]) > kCoefficientLimit) {
      return false;
    }
  }

  for (int sample = 0; sample < 720; ++sample) {
    const double theta =
        2.0 * kPi * static_cast<double>(sample) / 720.0;
    if (radial_factor(parameters, theta) < kMinimumRadialFactor) {
      return false;
    }
  }
  return true;
}

std::vector<Point2> shape_boundary(const ShapeParameters& parameters)
{
  if (!valid_shape(parameters)) {
    throw std::runtime_error("invalid Fourier shape parameters");
  }

  const double scale = area_scale(parameters);
  std::vector<Point2> points(kHoleBoundaryPointCount);
  for (int index = 0; index < kHoleBoundaryPointCount; ++index) {
    const double theta =
        2.0 * kPi * static_cast<double>(index)
        / static_cast<double>(kHoleBoundaryPointCount);
    const double radius = scale * radial_factor(parameters, theta);
    points[static_cast<std::size_t>(index)] = {
        radius * std::cos(theta),
        radius * std::sin(theta)};
  }

  // The Fourier scaling fixes the continuous polar area.  EasyMesh sees the
  // polygon through the sampled boundary points, so apply one final uniform
  // scaling that makes the discrete polygon area exactly kTargetHoleArea.
  double twice_polygon_area = 0.0;
  for (std::size_t index = 0; index < points.size(); ++index) {
    const Point2& point = points[index];
    const Point2& next = points[(index + 1) % points.size()];
    twice_polygon_area += point.x * next.y - next.x * point.y;
  }
  const double polygon_area = 0.5 * std::abs(twice_polygon_area);
  if (!(polygon_area > 0.0)) {
    throw std::runtime_error("degenerate sampled shape boundary");
  }
  const double polygon_scale =
      std::sqrt(kTargetHoleArea / polygon_area);
  for (Point2& point : points) {
    point.x *= polygon_scale;
    point.y *= polygon_scale;
  }
  return points;
}

std::vector<Point2> outer_boundary()
{
  std::vector<Point2> points(kOuterBoundaryPointCount);
  for (int index = 0; index < kOuterBoundaryPointCount; ++index) {
    const double theta =
        2.0 * kPi * static_cast<double>(index)
        / static_cast<double>(kOuterBoundaryPointCount);
    points[static_cast<std::size_t>(index)] = {
        kOuterRadius * std::cos(theta),
        kOuterRadius * std::sin(theta)};
  }
  return points;
}

void write_easymesh_domain(
    const fs::path& filename,
    const std::vector<Point2>& outer_points,
    const std::vector<Point2>& hole_points)
{
  std::ofstream output(filename);
  if (!output) {
    throw std::runtime_error(
        "cannot write EasyMesh domain " + filename.string());
  }

  output << std::setprecision(16);
  output << outer_points.size() + hole_points.size() << "\n";
  for (std::size_t index = 0; index < outer_points.size(); ++index) {
    const Point2& point = outer_points[index];
    output << index << ": "
           << point.x << " " << point.y << " "
           << kMeshSize << " " << kOuterBoundaryMark << "\n";
  }
  const std::size_t hole_offset = outer_points.size();
  for (std::size_t index = 0; index < hole_points.size(); ++index) {
    const Point2& point = hole_points[index];
    output << hole_offset + index << ": "
           << point.x << " " << point.y << " "
           << kMeshSize << " " << kHoleBoundaryMark << "\n";
  }

  output << outer_points.size() + hole_points.size() << "\n";
  for (std::size_t index = 0; index < outer_points.size(); ++index) {
    output << index << ": " << index << " "
           << (index + 1) % outer_points.size() << " "
           << kOuterBoundaryMark << "\n";
  }

  // A hole is represented by a clockwise boundary chain.
  for (std::size_t index = 0; index < hole_points.size(); ++index) {
    const std::size_t start =
        hole_offset + (hole_points.size() - index) % hole_points.size();
    const std::size_t end =
        hole_offset
        + (hole_points.size() - index - 1) % hole_points.size();
    output << outer_points.size() + index << ": "
           << start << " " << end << " "
           << kHoleBoundaryMark << "\n";
  }
}

void write_shape_csv(
    const fs::path& filename,
    const std::vector<Point2>& points)
{
  std::ofstream output(filename);
  output << "index,x,y\n";
  output << std::setprecision(16);
  for (std::size_t index = 0; index < points.size(); ++index) {
    output << index << "," << points[index].x << ","
           << points[index].y << "\n";
  }
}

void validate_easymesh_files(const fs::path& basename)
{
  const fs::path node_file = basename.string() + ".n";
  const fs::path element_file = basename.string() + ".e";
  const fs::path side_file = basename.string() + ".s";
  for (const fs::path& filename :
       {node_file, element_file, side_file}) {
    if (!fs::exists(filename) || fs::file_size(filename) == 0) {
      throw std::runtime_error(
          "EasyMesh did not create " + filename.string());
    }
  }

  std::ifstream nodes(node_file);
  std::ifstream elements(element_file);
  std::ifstream sides(side_file);
  int node_count = 0;
  int node_elements = 0;
  int node_sides = 0;
  int element_count = 0;
  int element_nodes = 0;
  int element_sides = 0;
  int side_count = 0;
  nodes >> node_count >> node_elements >> node_sides;
  elements >> element_count >> element_nodes >> element_sides;
  sides >> side_count;

  if (!nodes || !elements || !sides
      || node_count
          <= kOuterBoundaryPointCount + kHoleBoundaryPointCount
      || element_count <= 0
      || side_count <= 0
      || node_elements != element_count
      || node_sides != side_count
      || element_nodes != node_count
      || element_sides != side_count) {
    throw std::runtime_error(
        "inconsistent EasyMesh node/element/side headers");
  }
}

void run_easymesh(
    const fs::path& directory,
    const fs::path& easymesh_executable)
{
  const fs::path basename = directory / "shape";
  const std::string command =
      "cd " + shell_quote(directory.string())
      + " && " + shell_quote(easymesh_executable.string())
      + " shape.d > easymesh.log 2>&1";
  (void)std::system(command.c_str());
  validate_easymesh_files(basename);
}

void build_p1_space(
    EasyMesh& mesh,
    std::vector<TemplateElement<double, kDimension, kDimension> >&
        template_elements,
    FEMSpace<double, kDimension>& fem_space)
{
  fem_space.reinit(mesh, template_elements);
  const int element_count = mesh.n_geometry(kDimension);
  fem_space.element().resize(element_count);
  for (int index = 0; index < element_count; ++index) {
    if (mesh.geometry(kDimension, index).n_vertex() != 3) {
      throw std::runtime_error("the example expects triangle elements");
    }
    fem_space.element(index).reinit(fem_space, index, 0);
  }
  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();
}

void write_solution_svg(
    const fs::path& filename,
    EasyMesh& mesh,
    FEMSpace<double, kDimension>& fem_space,
    const FEMFunction<double, kDimension>& solution)
{
  struct TriangleValue {
    std::array<Point2, 3> vertices;
    double value = 0.0;
  };

  std::vector<TriangleValue> triangles;
  triangles.reserve(mesh.n_geometry(kDimension));
  double minimum_x = std::numeric_limits<double>::infinity();
  double maximum_x = -std::numeric_limits<double>::infinity();
  double minimum_y = std::numeric_limits<double>::infinity();
  double maximum_y = -std::numeric_limits<double>::infinity();
  double minimum_value = std::numeric_limits<double>::infinity();
  double maximum_value = -std::numeric_limits<double>::infinity();

  for (unsigned int index = 0;
       index < mesh.n_geometry(kDimension);
       ++index) {
    const Geometry& geometry = mesh.geometry(kDimension, index);
    TriangleValue triangle;
    Point<kDimension> centroid;
    centroid[0] = 0.0;
    centroid[1] = 0.0;
    for (int local = 0; local < 3; ++local) {
      const Point<kDimension>& point =
          mesh.point(geometry.vertex(local));
      triangle.vertices[static_cast<std::size_t>(local)] = {
          point[0], point[1]};
      centroid[0] += point[0] / 3.0;
      centroid[1] += point[1] / 3.0;
      minimum_x = std::min(minimum_x, point[0]);
      maximum_x = std::max(maximum_x, point[0]);
      minimum_y = std::min(minimum_y, point[1]);
      maximum_y = std::max(maximum_y, point[1]);
    }
    const std::vector<Point<kDimension> > query_points(1, centroid);
    triangle.value =
        solution.value(query_points, fem_space.element(index))[0];
    minimum_value = std::min(minimum_value, triangle.value);
    maximum_value = std::max(maximum_value, triangle.value);
    triangles.push_back(triangle);
  }

  constexpr double canvas_width = 760.0;
  constexpr double canvas_height = 620.0;
  constexpr double margin = 25.0;
  const double scale = std::min(
      (canvas_width - 2.0 * margin) / (maximum_x - minimum_x),
      (canvas_height - 2.0 * margin) / (maximum_y - minimum_y));
  const double offset_x =
      0.5 * (canvas_width - scale * (maximum_x - minimum_x));
  const double offset_y =
      0.5 * (canvas_height - scale * (maximum_y - minimum_y));
  const auto map_x = [&](double x) {
    return offset_x + scale * (x - minimum_x);
  };
  const auto map_y = [&](double y) {
    return canvas_height - offset_y - scale * (y - minimum_y);
  };

  std::ofstream output(filename);
  output << std::setprecision(8);
  output << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 " << canvas_width << " " << canvas_height
         << "\">\n";
  output << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
  const double range = std::max(1.0e-14, maximum_value - minimum_value);
  for (const TriangleValue& triangle : triangles) {
    const double normalized =
        (triangle.value - minimum_value) / range;
    const double hue = 240.0 * (1.0 - normalized);
    output << "<polygon points=\"";
    for (const Point2& point : triangle.vertices) {
      output << map_x(point.x) << "," << map_y(point.y) << " ";
    }
    output << "\" fill=\"hsl(" << hue
           << ",75%,55%)\" stroke=\"#425466\" "
           << "stroke-width=\"0.32\"/>\n";
  }
  output << "</svg>\n";
}

std::vector<Triangle2> extract_mesh_triangles(EasyMesh& mesh)
{
  std::vector<Triangle2> triangles;
  triangles.reserve(mesh.n_geometry(kDimension));
  for (unsigned int index = 0;
       index < mesh.n_geometry(kDimension);
       ++index) {
    const Geometry& geometry = mesh.geometry(kDimension, index);
    Triangle2 triangle;
    for (int local = 0; local < 3; ++local) {
      const Point<kDimension>& point =
          mesh.point(geometry.vertex(local));
      triangle[static_cast<std::size_t>(local)] = {
          point[0], point[1]};
    }
    triangles.push_back(triangle);
  }
  return triangles;
}

void write_mesh_svg(
    const fs::path& filename,
    const std::vector<Triangle2>& triangles)
{
  if (triangles.empty()) {
    return;
  }

  double minimum_x = std::numeric_limits<double>::infinity();
  double maximum_x = -std::numeric_limits<double>::infinity();
  double minimum_y = std::numeric_limits<double>::infinity();
  double maximum_y = -std::numeric_limits<double>::infinity();
  for (const Triangle2& triangle : triangles) {
    for (const Point2& point : triangle) {
      minimum_x = std::min(minimum_x, point.x);
      maximum_x = std::max(maximum_x, point.x);
      minimum_y = std::min(minimum_y, point.y);
      maximum_y = std::max(maximum_y, point.y);
    }
  }

  constexpr double width = 760.0;
  constexpr double height = 620.0;
  constexpr double margin = 25.0;
  const double scale = std::min(
      (width - 2.0 * margin) / (maximum_x - minimum_x),
      (height - 2.0 * margin) / (maximum_y - minimum_y));
  const double offset_x =
      0.5 * (width - scale * (maximum_x - minimum_x));
  const double offset_y =
      0.5 * (height - scale * (maximum_y - minimum_y));

  std::ofstream output(filename);
  output << std::setprecision(8)
         << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 " << width << " " << height << "\">\n"
         << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
  for (const Triangle2& triangle : triangles) {
    output << "<polygon points=\"";
    for (const Point2& point : triangle) {
      output << offset_x + scale * (point.x - minimum_x) << ","
             << height - offset_y - scale * (point.y - minimum_y)
             << " ";
    }
    output << "\" fill=\"#f8fafc\" stroke=\"#526b82\" "
           << "stroke-width=\"0.45\"/>\n";
  }
  output << "</svg>\n";
}

SolveResult solve_shape(
    const ShapeParameters& parameters,
    const fs::path& directory,
    const fs::path& easymesh_executable,
    std::vector<
        TemplateElement<double, kDimension, kDimension> >&
        template_elements,
    bool keep_artifacts)
{
  fs::remove_all(directory);
  fs::create_directories(directory);

  const std::vector<Point2> outer_points = outer_boundary();
  const std::vector<Point2> hole_points = shape_boundary(parameters);
  write_easymesh_domain(
      directory / "shape.d",
      outer_points,
      hole_points);
  write_shape_csv(directory / "outer_boundary.csv", outer_points);
  write_shape_csv(directory / "shape.csv", hole_points);
  run_easymesh(directory, easymesh_executable);

  EasyMesh mesh;
  const std::string basename = (directory / "shape").string();
  mesh.readData(basename.c_str());

  FEMSpace<double, kDimension> fem_space;
  build_p1_space(mesh, template_elements, fem_space);

  StiffMatrix<kDimension, double> stiffness(fem_space);
  stiffness.algebricAccuracy() = 4;
  stiffness.build();

  FEMFunction<double, kDimension> solution(fem_space);
  Vector<double> right_hand_side;
  Operator::L2Discretize(
      &unit_heat_source,
      fem_space,
      right_hand_side,
      4);

  BoundaryFunction<double, kDimension> cooled_outer_boundary(
      BoundaryConditionInfo::DIRICHLET,
      kOuterBoundaryMark,
      &zero_temperature);
  BoundaryFunction<double, kDimension> cooled_hole_boundary(
      BoundaryConditionInfo::DIRICHLET,
      kHoleBoundaryMark,
      &zero_temperature);
  BoundaryConditionAdmin<double, kDimension> boundary_admin(fem_space);
  boundary_admin.add(cooled_outer_boundary);
  boundary_admin.add(cooled_hole_boundary);
  boundary_admin.apply(
      stiffness,
      solution,
      right_hand_side);

  AMGSolver solver(stiffness);
  solver.solve(solution, right_hand_side, 1.0e-10, 300);

  SolveResult result;
  result.degrees_of_freedom = fem_space.n_dof();
  result.elements = mesh.n_geometry(kDimension);
  if (keep_artifacts) {
    result.mesh_triangles = extract_mesh_triangles(mesh);
  }
  for (unsigned int index = 0;
       index < mesh.n_geometry(kDimension);
       ++index) {
    const Element<double, kDimension>& element =
        fem_space.element(index);
    const QuadratureInfo<kDimension>& quadrature =
        element.findQuadratureInfo(4);
    const std::vector<Point<kDimension> > quadrature_points =
        element.local_to_global(quadrature.quadraturePoint());
    const std::vector<double> jacobians =
        element.local_to_global_jacobian(
            quadrature.quadraturePoint());
    const std::vector<double> values =
        solution.value(quadrature_points, element);
    const double reference_area = element.templateElement().volume();
    for (int q = 0; q < quadrature.n_quadraturePoint(); ++q) {
      const double weight =
          quadrature.weight(q) * jacobians[q] * reference_area;
      result.area += weight;
      result.temperature_integral += weight * values[q];
    }
  }
  result.mean_temperature =
      result.temperature_integral / result.area;

  if (keep_artifacts) {
    solution.writeOpenDXData((directory / "temperature.dx").string());
    write_mesh_svg(directory / "mesh.svg", result.mesh_triangles);
    write_solution_svg(
        directory / "temperature.svg",
        mesh,
        fem_space,
        solution);
    std::ofstream parameters_file(directory / "parameters.txt");
    parameters_file << std::setprecision(16)
                    << "a2 " << parameters.a2 << "\n"
                    << "b3 " << parameters.b3 << "\n"
                    << "a4 " << parameters.a4 << "\n"
                    << "domain_area " << result.area << "\n"
                    << "hole_area " << kTargetHoleArea << "\n"
                    << "mean_temperature "
                    << result.mean_temperature << "\n"
                    << "temperature_integral "
                    << result.temperature_integral << "\n"
                    << "dof " << result.degrees_of_freedom << "\n"
                    << "elements " << result.elements << "\n";
  }
  return result;
}

void write_outline_comparison(
    const fs::path& filename,
    const ShapeParameters& initial,
    const ShapeParameters& final_parameters)
{
  const std::vector<Point2> outer_points = outer_boundary();
  const std::vector<Point2> initial_points = shape_boundary(initial);
  const std::vector<Point2> final_points =
      shape_boundary(final_parameters);
  constexpr double width = 760.0;
  constexpr double height = 620.0;
  constexpr double scale = 135.0;
  constexpr double center_x = width / 2.0;
  constexpr double center_y = height / 2.0;

  const auto polygon_points =
      [&](const std::vector<Point2>& points) {
        std::ostringstream stream;
        stream << std::setprecision(8);
        for (const Point2& point : points) {
          stream << center_x + scale * point.x << ","
                 << center_y - scale * point.y << " ";
        }
        return stream.str();
      };

  std::ofstream output(filename);
  output << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 " << width << " " << height << "\">\n"
         << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n"
         << "<polygon points=\"" << polygon_points(outer_points)
         << "\" fill=\"#f8fafc\" stroke=\"#526b82\" "
         << "stroke-width=\"4\"/>\n"
         << "<polygon points=\"" << polygon_points(initial_points)
         << "\" fill=\"white\" stroke=\"#6b7280\" "
         << "stroke-width=\"4\" stroke-dasharray=\"10 8\"/>\n"
         << "<polygon points=\"" << polygon_points(final_points)
         << "\" fill=\"white\" stroke=\"#ea580c\" "
         << "stroke-width=\"5\"/>\n"
         << "</svg>\n";
}

void write_shape_history(
    const fs::path& filename,
    const std::vector<HistoryRecord>& records)
{
  if (records.empty()) {
    return;
  }

  constexpr double panel_width = 250.0;
  constexpr double height = 320.0;
  const double width = panel_width * records.size();
  const double center_y = 145.0;

  std::vector<std::vector<Point2> > boundaries;
  boundaries.reserve(records.size());
  const std::vector<Point2> outer_points = outer_boundary();
  const double maximum_extent = kOuterRadius;
  for (const HistoryRecord& record : records) {
    boundaries.push_back(shape_boundary(record.parameters));
  }

  const double scale = std::min(
      0.40 * panel_width / maximum_extent,
      0.36 * height / maximum_extent);

  std::ofstream output(filename);
  output << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 " << width << " " << height << "\">\n"
         << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

  for (std::size_t index = 0; index < records.size(); ++index) {
    const double center_x = panel_width * (index + 0.5);
    std::ostringstream outer_polygon;
    std::ostringstream hole_polygon;
    outer_polygon << std::setprecision(8);
    hole_polygon << std::setprecision(8);
    for (const Point2& point : outer_points) {
      outer_polygon << center_x + scale * point.x << ","
                    << center_y - scale * point.y << " ";
    }
    for (const Point2& point : boundaries[index]) {
      hole_polygon << center_x + scale * point.x << ","
                   << center_y - scale * point.y << " ";
    }

    output << "<polygon points=\"" << outer_polygon.str()
           << "\" fill=\"#f8fafc\" stroke=\"#526b82\" "
           << "stroke-width=\"2\"/>\n"
           << "<polygon points=\"" << hole_polygon.str()
           << "\" fill=\"white\" stroke=\"#ea580c\" "
           << "stroke-width=\"3\"/>\n"
           << "<text x=\"" << center_x
           << "\" y=\"275\" text-anchor=\"middle\" "
           << "font-family=\"Arial, sans-serif\" font-size=\"17\" "
           << "font-weight=\"600\" fill=\"#111827\">"
           << "iteration " << records[index].iteration
           << "</text>\n"
           << "<text x=\"" << center_x
           << "\" y=\"299\" text-anchor=\"middle\" "
           << "font-family=\"Arial, sans-serif\" font-size=\"15\" "
           << "fill=\"#4b5563\">J = "
           << std::fixed << std::setprecision(5)
           << records[index].result.mean_temperature
           << "</text>\n";
  }

  output << "</svg>\n";
}

void write_mesh_history(
    const fs::path& filename,
    const std::vector<HistoryRecord>& records)
{
  if (records.empty()) {
    return;
  }

  constexpr double panel_width = 250.0;
  constexpr double height = 320.0;
  const double width = panel_width * records.size();
  const double center_y = 145.0;

  double maximum_extent = 0.0;
  for (const HistoryRecord& record : records) {
    for (const Triangle2& triangle : record.result.mesh_triangles) {
      for (const Point2& point : triangle) {
        maximum_extent = std::max(
            maximum_extent,
            std::max(std::abs(point.x), std::abs(point.y)));
      }
    }
  }
  if (maximum_extent <= 0.0) {
    return;
  }

  const double scale = std::min(
      0.40 * panel_width / maximum_extent,
      0.36 * height / maximum_extent);

  std::ofstream output(filename);
  output << std::setprecision(8)
         << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         << "viewBox=\"0 0 " << width << " " << height << "\">\n"
         << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

  for (std::size_t index = 0; index < records.size(); ++index) {
    const double center_x = panel_width * (index + 0.5);
    for (const Triangle2& triangle
         : records[index].result.mesh_triangles) {
      output << "<polygon points=\"";
      for (const Point2& point : triangle) {
        output << center_x + scale * point.x << ","
               << center_y - scale * point.y << " ";
      }
      output << "\" fill=\"#f8fafc\" stroke=\"#526b82\" "
             << "stroke-width=\"0.32\"/>\n";
    }

    const std::vector<Point2> boundary =
        shape_boundary(records[index].parameters);
    output << "<polygon points=\"";
    for (const Point2& point : boundary) {
      output << center_x + scale * point.x << ","
             << center_y - scale * point.y << " ";
    }
    output << "\" fill=\"none\" stroke=\"#ea580c\" "
           << "stroke-width=\"1.8\"/>\n"
           << "<text x=\"" << center_x
           << "\" y=\"275\" text-anchor=\"middle\" "
           << "font-family=\"Arial, sans-serif\" font-size=\"17\" "
           << "font-weight=\"600\" fill=\"#111827\">"
           << "iteration " << records[index].iteration
           << "</text>\n"
           << "<text x=\"" << center_x
           << "\" y=\"299\" text-anchor=\"middle\" "
           << "font-family=\"Arial, sans-serif\" font-size=\"15\" "
           << "fill=\"#4b5563\">J = "
           << std::fixed << std::setprecision(5)
           << records[index].result.mean_temperature
           << "</text>\n";
  }

  output << "</svg>\n";
}

std::string iteration_name(int iteration)
{
  std::ostringstream name;
  name << "iter_" << std::setw(3) << std::setfill('0') << iteration;
  return name.str();
}

}  // namespace

int main(int argc, char** argv)
{
  try {
    if (argc > 3) {
      std::cerr << "usage: " << argv[0]
                << " [output_directory=output] [sweeps=4]\n";
      return 1;
    }

    const fs::path output_directory =
        argc > 1 ? fs::path(argv[1]) : fs::path("output");
    const int maximum_sweeps = argc > 2 ? std::stoi(argv[2]) : 4;
    if (maximum_sweeps < 1 || maximum_sweeps > 20) {
      throw std::runtime_error("sweeps must be in [1, 20]");
    }

    const fs::path home =
        environment_path("HOME", fs::current_path());
    const fs::path easymesh_executable =
        environment_path(
            "EASYMESH_BIN",
            home / "bin/easymesh");
    const fs::path template_directory =
        environment_path(
            "AFEPACK_TEMPLATE_PATH",
            home / "include/AFEPack/template/triangle");

    if (!fs::exists(easymesh_executable)) {
      throw std::runtime_error(
          "EasyMesh executable not found: "
          + easymesh_executable.string());
    }
    if (!fs::exists(template_directory / "triangle.tmp_geo")) {
      throw std::runtime_error(
          "AFEPack triangle templates not found: "
          + template_directory.string());
    }
    P1Templates templates;

    fs::create_directories(output_directory);
    fs::remove_all(output_directory / "iterations");
    fs::remove_all(output_directory / "work");
    fs::remove_all(output_directory / "reference_circle");
    fs::remove_all(output_directory / "reference_annulus");
    fs::create_directories(output_directory / "iterations");

    std::ofstream history_file(output_directory / "history.csv");
    history_file << "iteration,a2,b3,a4,step,domain_area,hole_area,"
                    "mean_temperature,temperature_integral,dof,elements\n";
    history_file << std::setprecision(16);
    std::vector<HistoryRecord> history_records;

    const ShapeParameters initial_parameters;
    ShapeParameters parameters = initial_parameters;
    double step = 0.05;

    SolveResult current = solve_shape(
        parameters,
        output_directory / "iterations" / iteration_name(0),
        easymesh_executable,
        templates.elements,
        true);
    history_file << 0 << "," << parameters.a2 << "," << parameters.b3
                 << "," << parameters.a4 << "," << step << ","
                 << current.area << "," << kTargetHoleArea << ","
                 << current.mean_temperature << ","
                 << current.temperature_integral << ","
                 << current.degrees_of_freedom << ","
                 << current.elements << "\n";
    history_records.push_back({0, parameters, current});

    std::cout << std::setprecision(10);
    std::cout << "Initial shape: a2=" << parameters.a2
              << ", b3=" << parameters.b3
              << ", a4=" << parameters.a4
              << ", mean T=" << current.mean_temperature << "\n";

    int completed_sweeps = 0;
    for (int sweep = 1; sweep <= maximum_sweeps; ++sweep) {
      bool improved_in_sweep = false;
      for (std::size_t coordinate = 0; coordinate < 3; ++coordinate) {
        ShapeParameters best_parameters = parameters;
        SolveResult best_result = current;

        for (double direction : {-1.0, 1.0}) {
          ShapeParameters trial = parameters;
          trial[coordinate] = std::clamp(
              trial[coordinate] + direction * step,
              -kCoefficientLimit,
              kCoefficientLimit);
          if (!valid_shape(trial)
              || trial[coordinate] == parameters[coordinate]) {
            continue;
          }

          try {
            const SolveResult trial_result = solve_shape(
                trial,
                output_directory / "work",
                easymesh_executable,
                templates.elements,
                false);
            std::cout << "  sweep " << sweep
                      << ", parameter " << coordinate
                      << ", trial=" << trial[coordinate]
                      << ", mean T="
                      << trial_result.mean_temperature << "\n";
            if (trial_result.mean_temperature
                > best_result.mean_temperature + 2.0e-7) {
              best_parameters = trial;
              best_result = trial_result;
            }
          } catch (const std::exception& error) {
            std::cout << "  rejected invalid trial: "
                      << error.what() << "\n";
          }
        }

        if (best_result.mean_temperature
            > current.mean_temperature + 2.0e-7) {
          parameters = best_parameters;
          current = best_result;
          improved_in_sweep = true;
        }
      }

      if (!improved_in_sweep) {
        step *= 0.5;
      }

      current = solve_shape(
          parameters,
          output_directory / "iterations" / iteration_name(sweep),
          easymesh_executable,
          templates.elements,
          true);
      history_file << sweep << "," << parameters.a2 << ","
                   << parameters.b3 << "," << parameters.a4 << ","
                   << step << "," << current.area << ","
                   << kTargetHoleArea << ","
                   << current.mean_temperature << ","
                   << current.temperature_integral << ","
                   << current.degrees_of_freedom << ","
                   << current.elements << "\n";
      history_file.flush();
      history_records.push_back({sweep, parameters, current});
      completed_sweeps = sweep;

      std::cout << "Accepted sweep " << sweep
                << ": a2=" << parameters.a2
                << ", b3=" << parameters.b3
                << ", a4=" << parameters.a4
                << ", step=" << step
                << ", mean T=" << current.mean_temperature << "\n";
      if (step < 0.0125) {
        break;
      }
    }

    ShapeParameters annulus_parameters;
    annulus_parameters.a2 = 0.0;
    annulus_parameters.b3 = 0.0;
    annulus_parameters.a4 = 0.0;
    const SolveResult annulus_result = solve_shape(
        annulus_parameters,
        output_directory / "reference_annulus",
        easymesh_executable,
        templates.elements,
        true);
    const double exact_annulus_mean =
        exact_concentric_annulus_mean_temperature();

    write_outline_comparison(
        output_directory / "shape_comparison.svg",
        initial_parameters,
        parameters);
    write_shape_history(
        output_directory / "shape_history.svg",
        history_records);
    write_mesh_history(
        output_directory / "mesh_history.svg",
        history_records);

    std::ofstream summary(output_directory / "summary.txt");
    summary << std::setprecision(16)
            << "completed_sweeps " << completed_sweeps << "\n"
            << "initial_a2 " << initial_parameters.a2 << "\n"
            << "initial_b3 " << initial_parameters.b3 << "\n"
            << "initial_a4 " << initial_parameters.a4 << "\n"
            << "final_a2 " << parameters.a2 << "\n"
            << "final_b3 " << parameters.b3 << "\n"
            << "final_a4 " << parameters.a4 << "\n"
            << "final_mean_temperature "
            << current.mean_temperature << "\n"
            << "discrete_concentric_annulus_mean_temperature "
            << annulus_result.mean_temperature << "\n"
            << "exact_concentric_annulus_mean_temperature "
            << exact_annulus_mean << "\n";

    fs::remove_all(output_directory / "work");

    std::cout << "\nFinal shape: a2=" << parameters.a2
              << ", b3=" << parameters.b3
              << ", a4=" << parameters.a4
              << ", mean T=" << current.mean_temperature << "\n";
    std::cout << "Discrete concentric annulus reference: mean T="
              << annulus_result.mean_temperature << "\n";
    std::cout << "Exact concentric annulus reference: mean T="
              << exact_annulus_mean << "\n";
    std::cout << "Outputs: " << fs::absolute(output_directory) << "\n";
  } catch (const std::exception& error) {
    std::cerr << "shape_optimization: " << error.what() << "\n";
    return 2;
  }
  return 0;
}
