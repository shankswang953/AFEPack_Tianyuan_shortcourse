#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <AFEPack/Geometry.h>
#include <AFEPack/HGeometry.h>

#include "local_mesh_merge.h"

namespace {

constexpr int kDimension = 2;
constexpr double kDomainWidth = 2.4;
constexpr double kDomainHeight = 1.0;

struct Stroke {
  double x0;
  double y0;
  double x1;
  double y1;
};

struct Glyph {
  const char* name;
  std::vector<Stroke> strokes;
};

Stroke shifted_stroke(double center_x,
                      double center_y,
                      double x0,
                      double y0,
                      double x1,
                      double y1) {
  return {
      center_x + x0,
      center_y + y0,
      center_x + x1,
      center_y + y1,
  };
}

Glyph make_tian_glyph() {
  constexpr double center_x = 0.65;
  constexpr double center_y = 0.50;
  return {
      "Tian",
      {
          // Tian: a short first horizontal, a longer second horizontal,
          // then a left-falling stroke that crosses the second horizontal
          // and a right-falling stroke beginning near that crossing.
          shifted_stroke(center_x, center_y, -0.25, 0.27, 0.25, 0.27),
          shifted_stroke(center_x, center_y, -0.43, 0.06, 0.43, 0.06),
          shifted_stroke(center_x, center_y, 0.05, 0.27, 0.01, 0.06),
          shifted_stroke(center_x, center_y, 0.01, 0.06, -0.38, -0.34),
          shifted_stroke(center_x, center_y, 0.01, 0.06, 0.40, -0.34),
      },
  };
}

Glyph make_yuan_glyph() {
  constexpr double center_x = 1.75;
  constexpr double center_y = 0.50;
  return {
      "Yuan",
      {
          shifted_stroke(center_x, center_y, -0.34, 0.25, 0.34, 0.25),
          shifted_stroke(center_x, center_y, -0.43, 0.06, 0.43, 0.06),
          shifted_stroke(center_x, center_y, -0.12, 0.06, -0.24, -0.33),
          shifted_stroke(center_x, center_y, 0.12, 0.06, 0.10, -0.26),
          shifted_stroke(center_x, center_y, 0.10, -0.26, 0.39, -0.26),
          shifted_stroke(center_x, center_y, 0.39, -0.26, 0.43, -0.16),
      },
  };
}

double distance_to_segment(double x, double y, const Stroke& stroke) {
  const double dx = stroke.x1 - stroke.x0;
  const double dy = stroke.y1 - stroke.y0;
  const double length_squared = dx * dx + dy * dy;
  if (length_squared == 0.0) {
    return std::hypot(x - stroke.x0, y - stroke.y0);
  }

  const double parameter = std::clamp(
      ((x - stroke.x0) * dx + (y - stroke.y0) * dy) / length_squared,
      0.0,
      1.0);
  const double closest_x = stroke.x0 + parameter * dx;
  const double closest_y = stroke.y0 + parameter * dy;
  return std::hypot(x - closest_x, y - closest_y);
}

// Distance to the closed thick-stroke set. It is exactly zero inside a stroke
// and positive outside, which gives the requested binary refinement rule.
double distance_to_glyph(double x,
                         double y,
                         const Glyph& glyph,
                         double half_width) {
  double centerline_distance = std::numeric_limits<double>::infinity();
  for (const Stroke& stroke : glyph.strokes) {
    centerline_distance =
        std::min(centerline_distance, distance_to_segment(x, y, stroke));
  }
  return std::max(centerline_distance - half_width, 0.0);
}

void element_barycenter(const RegularMesh<kDimension>& mesh,
                        int element_index,
                        double& x,
                        double& y) {
  const GeometryBM& element = mesh.geometry(kDimension, element_index);
  x = 0.0;
  y = 0.0;
  for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
    const Point<kDimension>& point = mesh.point(element.vertex(vertex));
    x += point[0];
    y += point[1];
  }
  x /= element.n_vertex();
  y /= element.n_vertex();
}

int set_glyph_indicator(RegularMesh<kDimension>& mesh,
                        Indicator<kDimension>& indicator,
                        const Glyph& glyph,
                        double half_width) {
  int selected = 0;
  for (int index = 0; index < mesh.n_geometry(kDimension); ++index) {
    indicator[index] = 0.0;
    double x = 0.0;
    double y = 0.0;
    element_barycenter(mesh, index, x, y);
    if (distance_to_glyph(x, y, glyph, half_width) == 0.0) {
      indicator[index] = 1.0;
      ++selected;
    }
  }
  return selected;
}

void refine_glyph(IrregularMesh<kDimension>& mesh,
                  const Glyph& glyph,
                  int rounds,
                  double half_width) {
  std::cout << "Refining the glyph " << glyph.name
            << " with d(x_K, S_" << glyph.name << ") = 0\n";

  for (int round = 0; round < rounds; ++round) {
    mesh.semiregularize();
    mesh.regularize(false);
    RegularMesh<kDimension>& regular_mesh = mesh.regularMesh();

    Indicator<kDimension> indicator(regular_mesh);
    const int selected =
        set_glyph_indicator(regular_mesh, indicator, glyph, half_width);
    if (selected == 0) {
      throw std::runtime_error(
          std::string("No element barycenter lies inside glyph ") + glyph.name);
    }

    std::cout << "  round " << round + 1
              << ": elements=" << regular_mesh.n_geometry(kDimension)
              << ", distance-zero elements=" << selected << '\n';

    MeshAdaptor<kDimension> adaptor(mesh);
    adaptor.convergenceOrder() = 0.0;
    adaptor.refineStep() = 0;
    adaptor.refineThreshold() = 1.0;
    adaptor.tolerence() = 1.0e-2;
    adaptor.is_refine_only() = true;
    adaptor.setIndicator(indicator);
    adaptor.adapt();
  }
}

void write_mesh(IrregularMesh<kDimension>& mesh,
                const std::string& basename) {
  mesh.semiregularize();
  mesh.regularize(false);
  RegularMesh<kDimension>& result = mesh.regularMesh();
  result.writeData(basename + ".mesh");
  result.writeEasyMesh(basename);
  result.writeOpenDXData(basename + ".dx");
  std::cout << "  " << basename << ": nodes=" << result.n_point()
            << ", regular elements=" << result.n_geometry(kDimension) << '\n';
}

void write_root_distance_field(HGeometryTree<kDimension>& hierarchy,
                               const Glyph& tian,
                               const Glyph& yuan,
                               double half_width) {
  IrregularMesh<kDimension> root_mesh(hierarchy);
  root_mesh.semiregularize();
  root_mesh.regularize(false);
  RegularMesh<kDimension>& mesh = root_mesh.regularMesh();
  mesh.writeOpenDXData("T_root.dx");

  std::ofstream output("glyph_distance_field.csv");
  output << std::setprecision(16);
  output << "element,x0,y0,x1,y1,x2,y2,center_x,center_y,"
            "distance_tian,distance_yuan,inside_tian,inside_yuan\n";

  for (int index = 0; index < mesh.n_geometry(kDimension); ++index) {
    const GeometryBM& element = mesh.geometry(kDimension, index);
    double x = 0.0;
    double y = 0.0;
    element_barycenter(mesh, index, x, y);
    const double distance_tian =
        distance_to_glyph(x, y, tian, half_width);
    const double distance_yuan =
        distance_to_glyph(x, y, yuan, half_width);

    output << index;
    for (int vertex = 0; vertex < 3; ++vertex) {
      const Point<kDimension>& point = mesh.point(element.vertex(vertex));
      output << ',' << point[0] << ',' << point[1];
    }
    output << ',' << x
           << ',' << y
           << ',' << distance_tian
           << ',' << distance_yuan
           << ',' << (distance_tian == 0.0 ? 1 : 0)
           << ',' << (distance_yuan == 0.0 ? 1 : 0)
           << '\n';
  }
}

void print_usage(const char* program) {
  std::cerr << "Usage: " << program
            << " [root_mesh_basename] [rounds] [stroke_half_width]\n"
            << "Default: " << program << " T_root 3 0.035\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc > 4) {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  const std::string root_basename = argc >= 2 ? argv[1] : "T_root";
  const int rounds = argc >= 3 ? std::atoi(argv[2]) : 3;
  const double half_width = argc >= 4 ? std::atof(argv[3]) : 0.035;
  if (rounds <= 0 || half_width <= 0.0 || half_width >= 0.12) {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  try {
    const Glyph tian = make_tian_glyph();
    const Glyph yuan = make_yuan_glyph();

    HGeometryTree<kDimension> hierarchy;
    hierarchy.readEasyMesh(root_basename);
    write_root_distance_field(hierarchy, tian, yuan, half_width);

    // Both peer meshes share exactly the same geometry tree. Their active
    // states differ, but root and child ordering remain compatible.
    IrregularMesh<kDimension> tian_mesh(hierarchy);
    IrregularMesh<kDimension> yuan_mesh(hierarchy);

    refine_glyph(tian_mesh, tian, rounds, half_width);
    std::cout << "\nTian input mesh:\n";
    write_mesh(tian_mesh, "T_tian");

    refine_glyph(yuan_mesh, yuan, rounds, half_width);
    std::cout << "\nYuan input mesh:\n";
    write_mesh(yuan_mesh, "T_yuan");

    local_multimesh::CommonIrregularMesh<kDimension> common_mesh(tian_mesh);
    common_mesh.merge(yuan_mesh);

    std::cout << "\nMerged Tianyuan mesh:\n";
    write_mesh(common_mesh, "T_common");
    std::cout
        << "\nThe common mesh contains both glyph refinement histories.\n"
        << "Domain: [0, " << kDomainWidth << "] x [0, "
        << kDomainHeight << "]\n";
  } catch (const std::exception& error) {
    std::cerr << "tianyuan_mesh_merge: " << error.what() << '\n';
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
