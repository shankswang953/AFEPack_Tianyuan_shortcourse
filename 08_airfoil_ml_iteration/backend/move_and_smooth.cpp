#include <AFEPack/HGeometry.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct Point2 {
  double x = 0.0;
  double y = 0.0;
};

struct Quality {
  double minimum_area = std::numeric_limits<double>::infinity();
  double minimum_shape_quality = std::numeric_limits<double>::infinity();
  int inverted_elements = 0;
};

std::vector<Point2> read_indexed_points(const fs::path& filename)
{
  std::ifstream input(filename);
  if (!input) {
    throw std::runtime_error("cannot open boundary points: " + filename.string());
  }

  std::vector<Point2> points;
  std::string line;
  while (std::getline(input, line)) {
    const std::size_t colon = line.find(':');
    if (colon == std::string::npos) {
      continue;
    }
    std::istringstream values(line.substr(colon + 1));
    Point2 point;
    if (!(values >> point.x >> point.y)) {
      throw std::runtime_error("invalid boundary point line");
    }
    points.push_back(point);
  }
  return points;
}

double twice_signed_area(
    const AFEPack::Point<2>& a,
    const AFEPack::Point<2>& b,
    const AFEPack::Point<2>& c)
{
  return (b[0] - a[0]) * (c[1] - a[1])
      - (b[1] - a[1]) * (c[0] - a[0]);
}

std::vector<double> element_orientations(const RegularMesh<2>& mesh)
{
  std::vector<double> orientation(
      static_cast<std::size_t>(mesh.n_geometry(2)), 0.0);
  for (int element = 0; element < mesh.n_geometry(2); ++element) {
    const Geometry& geometry = mesh.geometry(2, element);
    if (geometry.n_vertex() != 3) {
      throw std::runtime_error("the teaching example expects triangles");
    }
    const auto& vertices = geometry.vertex();
    orientation[static_cast<std::size_t>(element)] = twice_signed_area(
        mesh.point(vertices[0]),
        mesh.point(vertices[1]),
        mesh.point(vertices[2]));
  }
  return orientation;
}

Quality measure_quality(
    const RegularMesh<2>& mesh,
    const std::vector<double>& reference_orientation)
{
  Quality quality;
  for (int element = 0; element < mesh.n_geometry(2); ++element) {
    const Geometry& geometry = mesh.geometry(2, element);
    const auto& vertices = geometry.vertex();
    const AFEPack::Point<2>& a = mesh.point(vertices[0]);
    const AFEPack::Point<2>& b = mesh.point(vertices[1]);
    const AFEPack::Point<2>& c = mesh.point(vertices[2]);

    const double signed_twice_area = twice_signed_area(a, b, c);
    const double area = 0.5 * std::abs(signed_twice_area);
    const double ab2 = std::pow(a[0] - b[0], 2) + std::pow(a[1] - b[1], 2);
    const double bc2 = std::pow(b[0] - c[0], 2) + std::pow(b[1] - c[1], 2);
    const double ca2 = std::pow(c[0] - a[0], 2) + std::pow(c[1] - a[1], 2);
    const double shape_quality =
        4.0 * std::sqrt(3.0) * area / (ab2 + bc2 + ca2);

    quality.minimum_area = std::min(quality.minimum_area, area);
    quality.minimum_shape_quality =
        std::min(quality.minimum_shape_quality, shape_quality);
    if (signed_twice_area
        * reference_orientation[static_cast<std::size_t>(element)] <= 0.0) {
      ++quality.inverted_elements;
    }
  }
  return quality;
}

void write_mesh_csv(
    const RegularMesh<2>& mesh,
    const fs::path& output_dir,
    const std::string& stage)
{
  {
    std::ofstream nodes(output_dir / (stage + "_nodes.csv"));
    nodes << "index,x,y,boundary_mark\n";
    nodes << std::setprecision(16);
    for (int i = 0; i < mesh.n_point(); ++i) {
      nodes << i << "," << mesh.point(i)[0] << "," << mesh.point(i)[1]
            << "," << mesh.boundaryMark(0, i) << "\n";
    }
  }
  {
    std::ofstream elements(output_dir / (stage + "_elements.csv"));
    elements << "index,v0,v1,v2\n";
    for (int i = 0; i < mesh.n_geometry(2); ++i) {
      const Geometry& geometry = mesh.geometry(2, i);
      const auto& vertices = geometry.vertex();
      elements << i << "," << vertices[0] << "," << vertices[1] << ","
               << vertices[2] << "\n";
    }
  }
}

void write_stage(
    RegularMesh<2>& mesh,
    const fs::path& output_dir,
    const std::string& stage)
{
  mesh.writeOpenDXData((output_dir / (stage + ".dx")).string());
  mesh.writeData((output_dir / (stage + ".mesh")).string());
  write_mesh_csv(mesh, output_dir, stage);
}

int move_airfoil_boundary(
    RegularMesh<2>& mesh,
    const std::vector<Point2>& initial,
    const std::vector<Point2>& moved,
    double tolerance)
{
  if (initial.size() != moved.size()) {
    throw std::runtime_error("initial and moved boundary sizes differ");
  }

  std::vector<bool> vertex_used(
      static_cast<std::size_t>(mesh.n_point()), false);
  int moved_vertices = 0;
  for (std::size_t boundary_index = 0;
       boundary_index < initial.size();
       ++boundary_index) {
    int best_vertex = -1;
    double best_distance = std::numeric_limits<double>::infinity();
    for (int vertex = 0; vertex < mesh.n_point(); ++vertex) {
      if (vertex_used[static_cast<std::size_t>(vertex)]
          || mesh.boundaryMark(0, vertex) != 3) {
        continue;
      }
      const double current_distance = std::hypot(
          mesh.point(vertex)[0] - initial[boundary_index].x,
          mesh.point(vertex)[1] - initial[boundary_index].y);
      if (current_distance < best_distance) {
        best_distance = current_distance;
        best_vertex = vertex;
      }
    }

    if (best_vertex < 0 || best_distance > tolerance) {
      throw std::runtime_error(
          "could not match airfoil boundary point "
          + std::to_string(boundary_index)
          + "; nearest distance = "
          + std::to_string(best_distance));
    }

    mesh.point(best_vertex)[0] = moved[boundary_index].x;
    mesh.point(best_vertex)[1] = moved[boundary_index].y;
    vertex_used[static_cast<std::size_t>(best_vertex)] = true;
    ++moved_vertices;
  }
  return moved_vertices;
}

std::vector<std::vector<int>> build_vertex_adjacency(const RegularMesh<2>& mesh)
{
  std::vector<std::set<int>> adjacency_sets(
      static_cast<std::size_t>(mesh.n_point()));
  for (int element = 0; element < mesh.n_geometry(2); ++element) {
    const Geometry& geometry = mesh.geometry(2, element);
    const auto& vertices = geometry.vertex();
    for (int local_i = 0; local_i < geometry.n_vertex(); ++local_i) {
      for (int local_j = 0; local_j < geometry.n_vertex(); ++local_j) {
        if (local_i == local_j) {
          continue;
        }
        adjacency_sets[static_cast<std::size_t>(vertices[local_i])].insert(
            vertices[local_j]);
      }
    }
  }

  std::vector<std::vector<int>> adjacency(adjacency_sets.size());
  for (std::size_t i = 0; i < adjacency_sets.size(); ++i) {
    adjacency[i].assign(adjacency_sets[i].begin(), adjacency_sets[i].end());
  }
  return adjacency;
}

void laplacian_smooth(
    RegularMesh<2>& mesh,
    int iterations,
    double relaxation)
{
  if (iterations < 0 || !(relaxation > 0.0 && relaxation <= 1.0)) {
    throw std::runtime_error("invalid smoothing parameters");
  }

  const std::vector<std::vector<int>> adjacency =
      build_vertex_adjacency(mesh);
  std::vector<Point2> current(static_cast<std::size_t>(mesh.n_point()));
  std::vector<Point2> next = current;

  for (int iteration = 0; iteration < iterations; ++iteration) {
    for (int vertex = 0; vertex < mesh.n_point(); ++vertex) {
      current[static_cast<std::size_t>(vertex)] = {
          mesh.point(vertex)[0],
          mesh.point(vertex)[1]};
    }
    next = current;

    for (int vertex = 0; vertex < mesh.n_point(); ++vertex) {
      if (mesh.boundaryMark(0, vertex) != 0
          || adjacency[static_cast<std::size_t>(vertex)].empty()) {
        continue;
      }

      Point2 average;
      for (int neighbor : adjacency[static_cast<std::size_t>(vertex)]) {
        average.x += current[static_cast<std::size_t>(neighbor)].x;
        average.y += current[static_cast<std::size_t>(neighbor)].y;
      }
      const double count =
          static_cast<double>(adjacency[static_cast<std::size_t>(vertex)].size());
      average.x /= count;
      average.y /= count;

      next[static_cast<std::size_t>(vertex)].x =
          current[static_cast<std::size_t>(vertex)].x
          + relaxation * (average.x - current[static_cast<std::size_t>(vertex)].x);
      next[static_cast<std::size_t>(vertex)].y =
          current[static_cast<std::size_t>(vertex)].y
          + relaxation * (average.y - current[static_cast<std::size_t>(vertex)].y);
    }

    for (int vertex = 0; vertex < mesh.n_point(); ++vertex) {
      if (mesh.boundaryMark(0, vertex) == 0) {
        mesh.point(vertex)[0] = next[static_cast<std::size_t>(vertex)].x;
        mesh.point(vertex)[1] = next[static_cast<std::size_t>(vertex)].y;
      }
    }
  }
}

void print_quality(const std::string& stage, const Quality& quality)
{
  std::cout << std::left << std::setw(24) << stage
            << " min area = " << std::setw(12) << quality.minimum_area
            << " min quality = " << std::setw(12)
            << quality.minimum_shape_quality
            << " inverted = " << quality.inverted_elements << "\n";
}

int main(int argc, char** argv)
{
  try {
    if (argc < 5 || argc > 7) {
      std::cerr
          << "usage: move_and_smooth INITIAL_MESH INITIAL_BOUNDARY"
          << " MOVED_BOUNDARY OUTPUT_DIR"
          << " [iterations=50] [relaxation=0.45]\n";
      return 1;
    }

    const fs::path initial_mesh_file = argv[1];
    const fs::path initial_boundary_file = argv[2];
    const fs::path moved_boundary_file = argv[3];
    const fs::path output_dir = argv[4];
    const int iterations = argc > 5 ? std::stoi(argv[5]) : 50;
    const double relaxation = argc > 6 ? std::stod(argv[6]) : 0.45;

    fs::create_directories(output_dir);

    HGeometryTree<2> geometry_tree;
    geometry_tree.readMesh(initial_mesh_file.string());
    IrregularMesh<2> irregular_mesh(geometry_tree);
    irregular_mesh.semiregularize();
    irregular_mesh.regularize(false);
    RegularMesh<2>& mesh = irregular_mesh.regularMesh();

    const std::vector<double> reference_orientation =
        element_orientations(mesh);
    const Quality initial_quality =
        measure_quality(mesh, reference_orientation);
    write_stage(mesh, output_dir, "mesh_initial");

    const std::vector<Point2> initial_boundary =
        read_indexed_points(initial_boundary_file);
    const std::vector<Point2> moved_boundary =
        read_indexed_points(moved_boundary_file);
    const int moved_vertices = move_airfoil_boundary(
        mesh,
        initial_boundary,
        moved_boundary,
        1.0e-7);

    const Quality moved_quality =
        measure_quality(mesh, reference_orientation);
    write_stage(mesh, output_dir, "mesh_moved_unsmoothed");

    laplacian_smooth(mesh, iterations, relaxation);
    const Quality smoothed_quality =
        measure_quality(mesh, reference_orientation);
    write_stage(mesh, output_dir, "mesh_smoothed");

    std::ofstream summary(output_dir / "quality_summary.csv");
    summary << "stage,min_area,min_shape_quality,inverted_elements\n";
    summary << std::setprecision(16);
    summary << "initial," << initial_quality.minimum_area << ","
            << initial_quality.minimum_shape_quality << ","
            << initial_quality.inverted_elements << "\n";
    summary << "moved_unsmoothed," << moved_quality.minimum_area << ","
            << moved_quality.minimum_shape_quality << ","
            << moved_quality.inverted_elements << "\n";
    summary << "smoothed," << smoothed_quality.minimum_area << ","
            << smoothed_quality.minimum_shape_quality << ","
            << smoothed_quality.inverted_elements << "\n";

    std::cout << "Mesh points: " << mesh.n_point()
              << ", triangle elements: " << mesh.n_geometry(2) << "\n";
    std::cout << "Matched and moved airfoil vertices: "
              << moved_vertices << "\n";
    print_quality("initial", initial_quality);
    print_quality("moved, unsmoothed", moved_quality);
    print_quality("after Laplacian smooth", smoothed_quality);
    if (smoothed_quality.inverted_elements != 0
        || !(smoothed_quality.minimum_area > 0.0)) {
      throw std::runtime_error(
          "smoothed mesh is invalid; persistent mesh was not updated");
    }
  } catch (const std::exception& error) {
    std::cerr << "move_and_smooth: " << error.what() << "\n";
    return 2;
  }
  return 0;
}
