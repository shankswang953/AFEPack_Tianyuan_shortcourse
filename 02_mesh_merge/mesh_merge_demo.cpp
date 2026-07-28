#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>

#include <AFEPack/Geometry.h>
#include <AFEPack/HGeometry.h>

#include "local_mesh_merge.h"

namespace {

constexpr int kDimension = 2;

struct RefinementRegion {
  const char* name;
  double center_x;
  double center_y;
  double radius;
};

int set_circular_indicator(RegularMesh<kDimension>& mesh,
                           Indicator<kDimension>& indicator,
                           const RefinementRegion& region) {
  int selected = 0;
  for (int i = 0; i < mesh.n_geometry(kDimension); ++i) {
    const GeometryBM& element = mesh.geometry(kDimension, i);
    double x = 0.0;
    double y = 0.0;
    for (int j = 0; j < element.n_vertex(); ++j) {
      const Point<kDimension>& point = mesh.point(element.vertex(j));
      x += point[0];
      y += point[1];
    }
    x /= element.n_vertex();
    y /= element.n_vertex();

    const double dx = x - region.center_x;
    const double dy = y - region.center_y;
    if (dx * dx + dy * dy <= region.radius * region.radius) {
      indicator[i] = 1.0;
      ++selected;
    }
  }
  return selected;
}

void refine_region(IrregularMesh<kDimension>& mesh,
                   const RefinementRegion& region,
                   int rounds) {
  std::cout << "Refining " << region.name << " around ("
            << region.center_x << ", " << region.center_y << ")\n";

  for (int round = 0; round < rounds; ++round) {
    mesh.semiregularize();
    mesh.regularize(false);
    RegularMesh<kDimension>& regular_mesh = mesh.regularMesh();

    Indicator<kDimension> indicator(regular_mesh);
    const int selected =
        set_circular_indicator(regular_mesh, indicator, region);
    if (selected == 0) {
      throw std::runtime_error(
          std::string("No element selected in region ") + region.name);
    }

    std::cout << "  round " << round + 1
              << ": elements=" << regular_mesh.n_geometry(kDimension)
              << ", selected=" << selected << '\n';

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
            << ", elements=" << result.n_geometry(kDimension) << '\n';
}

void print_usage(const char* program) {
  std::cerr << "Usage: " << program << " [root_mesh_basename] [rounds]\n"
            << "Default: " << program << " D 4\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc > 3) {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  const std::string root_basename = argc >= 2 ? argv[1] : "D";
  const int rounds = argc == 3 ? std::atoi(argv[2]) : 4;
  if (rounds <= 0) {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  const RefinementRegion lower_left = {
      "lower-left region", 0.18, 0.20, 0.11};
  const RefinementRegion upper_right = {
      "upper-right region", 0.82, 0.76, 0.11};

  try {
    // The two independent refinement states deliberately share one HGT.
    // This guarantees identical root and child ordering for the merge.
    HGeometryTree<kDimension> hierarchy;
    hierarchy.readEasyMesh(root_basename);
    IrregularMesh<kDimension> left_mesh(hierarchy);
    IrregularMesh<kDimension> right_mesh(hierarchy);

    refine_region(left_mesh, lower_left, rounds);
    std::cout << "\nLeft input mesh:\n";
    write_mesh(left_mesh, "D_left");

    refine_region(right_mesh, upper_right, rounds);
    std::cout << "\nRight input mesh:\n";
    write_mesh(right_mesh, "D_right");

    local_multimesh::CommonIrregularMesh<kDimension> common_mesh(left_mesh);
    common_mesh.merge(right_mesh);

    std::cout << "\nMerged common mesh:\n";
    write_mesh(common_mesh, "D_common");
    std::cout << "\nThe common mesh contains the union of both refinement trees.\n";
  } catch (const std::exception& error) {
    std::cerr << "mesh_merge_demo: " << error.what() << '\n';
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
