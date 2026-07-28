#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include <AFEPack/Geometry.h>
#include <AFEPack/HGeometry.h>

namespace {

constexpr int kDimension = 2;

int set_circular_indicator(RegularMesh<kDimension>& mesh,
                           Indicator<kDimension>& indicator,
                           double center_x,
                           double center_y,
                           double radius) {
  int selected = 0;

  for (int i = 0; i < mesh.n_geometry(kDimension); ++i) {
    const GeometryBM& element = mesh.geometry(kDimension, i);
    const int n_vertex = element.n_vertex();

    double barycenter_x = 0.0;
    double barycenter_y = 0.0;
    for (int j = 0; j < n_vertex; ++j) {
      const Point<kDimension>& point = mesh.point(element.vertex(j));
      barycenter_x += point[0];
      barycenter_y += point[1];
    }
    barycenter_x /= n_vertex;
    barycenter_y /= n_vertex;

    const double dx = barycenter_x - center_x;
    const double dy = barycenter_y - center_y;
    if (std::sqrt(dx * dx + dy * dy) <= radius) {
      indicator[i] = 1.0;
      ++selected;
    }
  }

  return selected;
}

void print_usage(const char* program) {
  std::cerr << "Usage: " << program
            << " input_basename output_basename center_x center_y radius rounds\n"
            << "Example: " << program
            << " D D_refine_lower_left 0.18 0.20 0.11 4\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 7) {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  const std::string input_basename = argv[1];
  const std::string output_basename = argv[2];
  const double center_x = std::atof(argv[3]);
  const double center_y = std::atof(argv[4]);
  const double radius = std::atof(argv[5]);
  const int rounds = std::atoi(argv[6]);

  if (radius <= 0.0 || rounds <= 0) {
    std::cerr << "radius and rounds must both be positive.\n";
    return EXIT_FAILURE;
  }

  HGeometryTree<kDimension> hierarchy;
  hierarchy.readEasyMesh(input_basename);
  IrregularMesh<kDimension> irregular_mesh(hierarchy);

  std::cout << "Input root mesh: " << input_basename << '\n'
            << "Indicator center: (" << center_x << ", " << center_y << ")\n"
            << "Indicator radius: " << radius << '\n'
            << "Refinement rounds: " << rounds << "\n\n";

  for (int round = 0; round < rounds; ++round) {
    irregular_mesh.semiregularize();
    irregular_mesh.regularize(false);
    RegularMesh<kDimension>& mesh = irregular_mesh.regularMesh();

    Indicator<kDimension> indicator(mesh);
    const int selected = set_circular_indicator(
        mesh, indicator, center_x, center_y, radius);

    std::cout << "Round " << round + 1
              << ": nodes=" << mesh.n_point()
              << ", elements=" << mesh.n_geometry(kDimension)
              << ", indicator-selected=" << selected << '\n';

    if (selected == 0) {
      std::cerr << "No element barycenter lies inside the indicator region.\n";
      return EXIT_FAILURE;
    }

    MeshAdaptor<kDimension> adaptor(irregular_mesh);
    adaptor.convergenceOrder() = 0.0;
    adaptor.refineStep() = 0;  // Exactly one local level per outer round.
    adaptor.refineThreshold() = 1.0;
    adaptor.tolerence() = 1.0e-2;
    adaptor.is_refine_only() = true;
    adaptor.setIndicator(indicator);
    adaptor.adapt();
  }

  irregular_mesh.semiregularize();
  irregular_mesh.regularize(false);
  RegularMesh<kDimension>& result = irregular_mesh.regularMesh();

  result.writeData(output_basename + ".mesh");
  result.writeEasyMesh(output_basename);
  result.writeOpenDXData(output_basename + ".dx");

  std::cout << "\nFinal mesh: nodes=" << result.n_point()
            << ", elements=" << result.n_geometry(kDimension) << '\n'
            << "Wrote " << output_basename << ".mesh, "
            << output_basename << ".[nse], and "
            << output_basename << ".dx\n";

  return EXIT_SUCCESS;
}
