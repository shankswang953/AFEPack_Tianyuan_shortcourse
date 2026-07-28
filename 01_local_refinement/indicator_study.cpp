#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <AFEPack/Geometry.h>
#include <AFEPack/HGeometry.h>

namespace {

constexpr int kDimension = 2;
constexpr double kCenterAx = 0.18;
constexpr double kCenterAy = 0.20;
constexpr double kCenterBx = 0.82;
constexpr double kCenterBy = 0.76;
constexpr double kSigma = 0.10;

struct ElementScore {
  int index = -1;
  double center_x = 0.0;
  double center_y = 0.0;
  double eta_a = 0.0;
  double eta_b = 0.0;
  double score = 0.0;
};

enum class IndicatorMode {
  kIndependentA,
  kIndependentB,
  kEqualSum,
  kScaledSum,
};

struct StudyCase {
  const char* key;
  const char* output_name;
  IndicatorMode mode;
  double marking_fraction;
};

double gaussian_indicator(double x,
                          double y,
                          double center_x,
                          double center_y) {
  const double dx = x - center_x;
  const double dy = y - center_y;
  return std::exp(-(dx * dx + dy * dy) / (2.0 * kSigma * kSigma));
}

ElementScore evaluate_element(const RegularMesh<kDimension>& mesh,
                              int element_index,
                              IndicatorMode mode,
                              double scale_b) {
  const GeometryBM& element = mesh.geometry(kDimension, element_index);
  ElementScore result;
  result.index = element_index;

  for (int vertex = 0; vertex < element.n_vertex(); ++vertex) {
    const Point<kDimension>& point = mesh.point(element.vertex(vertex));
    result.center_x += point[0];
    result.center_y += point[1];
  }
  result.center_x /= element.n_vertex();
  result.center_y /= element.n_vertex();

  result.eta_a = gaussian_indicator(
      result.center_x, result.center_y, kCenterAx, kCenterAy);
  result.eta_b = gaussian_indicator(
      result.center_x, result.center_y, kCenterBx, kCenterBy);

  switch (mode) {
    case IndicatorMode::kIndependentA:
      result.score = result.eta_a;
      break;
    case IndicatorMode::kIndependentB:
      result.score = result.eta_b;
      break;
    case IndicatorMode::kEqualSum:
      result.score = result.eta_a + result.eta_b;
      break;
    case IndicatorMode::kScaledSum:
      result.score = result.eta_a + scale_b * result.eta_b;
      break;
  }
  return result;
}

void write_root_indicator_fields(const std::string& input_basename,
                                 const std::filesystem::path& output_dir,
                                 double scale_b) {
  HGeometryTree<kDimension> hierarchy;
  hierarchy.readEasyMesh(input_basename);
  IrregularMesh<kDimension> irregular_mesh(hierarchy);
  irregular_mesh.semiregularize();
  irregular_mesh.regularize(false);
  RegularMesh<kDimension>& mesh = irregular_mesh.regularMesh();

  const std::string root_name = (output_dir / "D_root").string();
  mesh.writeData(root_name + ".mesh");
  mesh.writeEasyMesh(root_name);
  mesh.writeOpenDXData(root_name + ".dx");

  std::ofstream output(output_dir / "indicator_fields.csv");
  output << std::setprecision(16);
  output << "element,x0,y0,x1,y1,x2,y2,center_x,center_y,"
            "eta_a,eta_b,eta_equal,eta_scaled\n";

  for (int index = 0; index < mesh.n_geometry(kDimension); ++index) {
    const GeometryBM& element = mesh.geometry(kDimension, index);
    const ElementScore values = evaluate_element(
        mesh, index, IndicatorMode::kEqualSum, scale_b);
    output << index;
    for (int vertex = 0; vertex < 3; ++vertex) {
      const Point<kDimension>& point = mesh.point(element.vertex(vertex));
      output << ',' << point[0] << ',' << point[1];
    }
    output << ',' << values.center_x
           << ',' << values.center_y
           << ',' << values.eta_a
           << ',' << values.eta_b
           << ',' << values.eta_a + values.eta_b
           << ',' << values.eta_a + scale_b * values.eta_b
           << '\n';
  }
}

void run_case(const std::string& input_basename,
              const std::filesystem::path& output_dir,
              const StudyCase& study_case,
              int rounds,
              double scale_b,
              std::ofstream& history) {
  HGeometryTree<kDimension> hierarchy;
  hierarchy.readEasyMesh(input_basename);
  IrregularMesh<kDimension> irregular_mesh(hierarchy);

  std::cout << "\n[" << study_case.key << "] marking fraction = "
            << study_case.marking_fraction << '\n';

  for (int round = 0; round < rounds; ++round) {
    irregular_mesh.semiregularize();
    irregular_mesh.regularize(false);
    RegularMesh<kDimension>& mesh = irregular_mesh.regularMesh();
    const int element_count = mesh.n_geometry(kDimension);
    const int marked_count = std::max(
        1,
        std::min(
            element_count,
            static_cast<int>(
                std::ceil(study_case.marking_fraction * element_count))));

    std::vector<ElementScore> scores;
    scores.reserve(element_count);
    for (int index = 0; index < element_count; ++index) {
      scores.push_back(evaluate_element(
          mesh, index, study_case.mode, scale_b));
    }
    std::partial_sort(
        scores.begin(),
        scores.begin() + marked_count,
        scores.end(),
        [](const ElementScore& left, const ElementScore& right) {
          return left.score > right.score;
        });

    Indicator<kDimension> indicator(mesh);
    for (int index = 0; index < element_count; ++index) {
      indicator[index] = 0.0;
    }

    int marked_near_a = 0;
    int marked_near_b = 0;
    for (int rank = 0; rank < marked_count; ++rank) {
      const ElementScore& selected = scores[rank];
      indicator[selected.index] = 1.0;
      if (selected.eta_a >= selected.eta_b) {
        ++marked_near_a;
      } else {
        ++marked_near_b;
      }
    }

    const double cutoff = scores[marked_count - 1].score;
    std::cout << "  round " << round + 1
              << ": elements=" << element_count
              << ", marked=" << marked_count
              << " (A=" << marked_near_a
              << ", B=" << marked_near_b << ")"
              << ", cutoff=" << cutoff << '\n';
    history << study_case.key << ','
            << round + 1 << ','
            << element_count << ','
            << marked_count << ','
            << marked_near_a << ','
            << marked_near_b << ','
            << cutoff << ','
            << scores.front().score << '\n';

    MeshAdaptor<kDimension> adaptor(irregular_mesh);
    adaptor.convergenceOrder() = 0.0;
    adaptor.refineStep() = 0;
    adaptor.refineThreshold() = 1.0;
    adaptor.tolerence() = 1.0e-2;
    adaptor.is_refine_only() = true;
    adaptor.setIndicator(indicator);
    adaptor.adapt();
  }

  irregular_mesh.semiregularize();
  irregular_mesh.regularize(false);
  RegularMesh<kDimension>& result = irregular_mesh.regularMesh();
  const std::string output_name =
      (output_dir / study_case.output_name).string();
  result.writeData(output_name + ".mesh");
  result.writeEasyMesh(output_name);
  result.writeOpenDXData(output_name + ".dx");

  std::cout << "  final: nodes=" << result.n_point()
            << ", elements=" << result.n_geometry(kDimension)
            << ", output=" << study_case.output_name << '\n';
}

void print_usage(const char* program) {
  std::cerr
      << "Usage: " << program
      << " ROOT_MESH OUTPUT_DIR [ROUNDS] [LOCAL_FRACTION] [SCALE_B]\n"
      << "Defaults: rounds=3, local_fraction=0.05, scale_b=50\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc < 3 || argc > 6) {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  const std::string root_mesh = argv[1];
  const std::filesystem::path output_dir = argv[2];
  const int rounds = argc >= 4 ? std::atoi(argv[3]) : 3;
  const double local_fraction = argc >= 5 ? std::atof(argv[4]) : 0.05;
  const double scale_b = argc >= 6 ? std::atof(argv[5]) : 50.0;

  if (rounds <= 0 ||
      local_fraction <= 0.0 ||
      local_fraction >= 0.5 ||
      scale_b <= 0.0) {
    std::cerr
        << "Require rounds > 0, 0 < local_fraction < 0.5, and scale_b > 0.\n";
    return EXIT_FAILURE;
  }

  std::filesystem::create_directories(output_dir);
  write_root_indicator_fields(root_mesh, output_dir, scale_b);

  std::ofstream history(output_dir / "marking_history.csv");
  history << std::setprecision(16);
  history << "case,round,elements_before,marked_total,"
             "marked_near_a,marked_near_b,cutoff,max_score\n";

  const std::vector<StudyCase> cases = {
      {"independent_A",
       "D_refine_lower_left",
       IndicatorMode::kIndependentA,
       local_fraction},
      {"independent_B",
       "D_refine_upper_right",
       IndicatorMode::kIndependentB,
       local_fraction},
      {"equal_sum",
       "D_combined_equal",
       IndicatorMode::kEqualSum,
       2.0 * local_fraction},
      {"scaled_sum",
       "D_combined_scaled",
       IndicatorMode::kScaledSum,
       2.0 * local_fraction},
  };

  std::cout << std::setprecision(6)
            << "Two-indicator refinement study\n"
            << "  A center          : (" << kCenterAx << ", "
            << kCenterAy << ")\n"
            << "  B center          : (" << kCenterBx << ", "
            << kCenterBy << ")\n"
            << "  Gaussian sigma    : " << kSigma << '\n'
            << "  rounds            : " << rounds << '\n'
            << "  local fraction    : " << local_fraction << '\n'
            << "  combined fraction : " << 2.0 * local_fraction << '\n'
            << "  B scale coefficient: " << scale_b << '\n';

  for (const StudyCase& study_case : cases) {
    run_case(
        root_mesh,
        output_dir,
        study_case,
        rounds,
        scale_b,
        history);
  }

  std::cout << "\nStudy outputs are in " << output_dir << '\n';
  return EXIT_SUCCESS;
}
