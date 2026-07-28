#ifndef AIRFOIL_BEZIER_HPP
#define AIRFOIL_BEZIER_HPP

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <lapacke.h>

namespace airfoil_demo {

// This header is the teaching-sized extraction of the airfoil geometry part
// of RL/include/circleAirfoil.h.  It keeps the regularized Bezier fitting and
// boundary resampling ideas, but deliberately contains no Euler solver,
// optimization agent, or mesh-adaptation code.
//
// Legacy-name map:
//   LSFittingSmooth()        -> fit_regularized_bezier()
//   setup_piecewiseCurve()   -> BezierCurve
//   generate_easymeshPoints()-> fit_and_sample()
//
// The legacy file solves its least-squares system with GSL.  This standalone
// example uses LAPACKE_dgels so it builds with the numerical libraries already
// used by the local AFEPack installation.

struct Point2 {
  double x = 0.0;
  double y = 0.0;
};

struct AirfoilData {
  std::string title;
  std::vector<Point2> upper;
  std::vector<Point2> lower;
};

inline AirfoilData read_uiuc_airfoil(const std::filesystem::path& filename)
{
  std::ifstream input(filename);
  if (!input) {
    throw std::runtime_error("cannot open airfoil data: " + filename.string());
  }

  AirfoilData data;
  std::getline(input, data.title);

  std::string count_line;
  std::getline(input, count_line);
  std::istringstream count_stream(count_line);
  double upper_count_value = 0.0;
  double lower_count_value = 0.0;
  if (!(count_stream >> upper_count_value >> lower_count_value)) {
    throw std::runtime_error("invalid point-count line in " + filename.string());
  }

  const int upper_count = static_cast<int>(std::lround(upper_count_value));
  const int lower_count = static_cast<int>(std::lround(lower_count_value));
  if (upper_count < 2 || lower_count < 2) {
    throw std::runtime_error("an airfoil surface needs at least two points");
  }

  data.upper.resize(static_cast<std::size_t>(upper_count));
  data.lower.resize(static_cast<std::size_t>(lower_count));
  for (Point2& point : data.upper) {
    if (!(input >> point.x >> point.y)) {
      throw std::runtime_error("not enough upper-surface points");
    }
  }
  for (Point2& point : data.lower) {
    if (!(input >> point.x >> point.y)) {
      throw std::runtime_error("not enough lower-surface points");
    }
  }
  return data;
}

inline void write_uiuc_airfoil(
    const std::filesystem::path& filename,
    const AirfoilData& data)
{
  std::ofstream output(filename);
  if (!output) {
    throw std::runtime_error("cannot write airfoil data: " + filename.string());
  }

  output << data.title << "\n";
  output << data.upper.size() << " " << data.lower.size() << ".\n\n";
  output << std::setprecision(16) << std::fixed;
  for (const Point2& point : data.upper) {
    output << point.x << " " << point.y << "\n";
  }
  output << "\n";
  for (const Point2& point : data.lower) {
    output << point.x << " " << point.y << "\n";
  }
}

inline void close_airfoil_endpoints(AirfoilData& data)
{
  if (data.upper.empty() || data.lower.empty()) {
    throw std::runtime_error("cannot close an empty airfoil");
  }

  const Point2 leading = {
      0.5 * (data.upper.front().x + data.lower.front().x),
      0.5 * (data.upper.front().y + data.lower.front().y)};
  const Point2 trailing = {
      0.5 * (data.upper.back().x + data.lower.back().x),
      0.5 * (data.upper.back().y + data.lower.back().y)};
  data.upper.front() = leading;
  data.lower.front() = leading;
  data.upper.back() = trailing;
  data.lower.back() = trailing;
}

inline double binomial_coefficient(int n, int k)
{
  if (k < 0 || k > n) {
    return 0.0;
  }
  k = std::min(k, n - k);
  double result = 1.0;
  for (int i = 1; i <= k; ++i) {
    result *= static_cast<double>(n - k + i) / static_cast<double>(i);
  }
  return result;
}

inline double bernstein(int index, int degree, double t)
{
  if (t <= 0.0) {
    return index == 0 ? 1.0 : 0.0;
  }
  if (t >= 1.0) {
    return index == degree ? 1.0 : 0.0;
  }
  return binomial_coefficient(degree, index)
      * std::pow(t, index)
      * std::pow(1.0 - t, degree - index);
}

class BezierCurve {
 public:
  explicit BezierCurve(std::vector<Point2> control_points)
      : control_points_(std::move(control_points))
  {
    if (control_points_.size() < 2) {
      throw std::runtime_error("a Bezier curve needs at least two control points");
    }
  }

  const std::vector<Point2>& control_points() const
  {
    return control_points_;
  }

  Point2 evaluate(double t) const
  {
    std::vector<Point2> work = control_points_;
    for (std::size_t level = work.size() - 1; level > 0; --level) {
      for (std::size_t i = 0; i < level; ++i) {
        work[i].x = (1.0 - t) * work[i].x + t * work[i + 1].x;
        work[i].y = (1.0 - t) * work[i].y + t * work[i + 1].y;
      }
    }
    return work.front();
  }

  std::vector<Point2> sample_uniform_arclength(
      int count,
      int dense_count = 8001) const
  {
    if (count < 2 || dense_count < count) {
      throw std::runtime_error("invalid Bezier sampling size");
    }

    std::vector<Point2> dense(static_cast<std::size_t>(dense_count));
    std::vector<double> cumulative(static_cast<std::size_t>(dense_count), 0.0);
    for (int i = 0; i < dense_count; ++i) {
      const double t = static_cast<double>(i) / static_cast<double>(dense_count - 1);
      dense[static_cast<std::size_t>(i)] = evaluate(t);
      if (i > 0) {
        const Point2& p0 = dense[static_cast<std::size_t>(i - 1)];
        const Point2& p1 = dense[static_cast<std::size_t>(i)];
        cumulative[static_cast<std::size_t>(i)] =
            cumulative[static_cast<std::size_t>(i - 1)]
            + std::hypot(p1.x - p0.x, p1.y - p0.y);
      }
    }

    const double total_length = cumulative.back();
    if (!(total_length > 0.0)) {
      throw std::runtime_error("zero-length Bezier curve");
    }

    std::vector<Point2> result;
    result.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
      const double target =
          total_length * static_cast<double>(i) / static_cast<double>(count - 1);
      const auto upper = std::lower_bound(cumulative.begin(), cumulative.end(), target);
      if (upper == cumulative.begin()) {
        result.push_back(dense.front());
        continue;
      }
      if (upper == cumulative.end()) {
        result.push_back(dense.back());
        continue;
      }

      const std::size_t right = static_cast<std::size_t>(
          std::distance(cumulative.begin(), upper));
      const std::size_t left = right - 1;
      const double span = cumulative[right] - cumulative[left];
      const double alpha = span > 0.0 ? (target - cumulative[left]) / span : 0.0;
      result.push_back({
          (1.0 - alpha) * dense[left].x + alpha * dense[right].x,
          (1.0 - alpha) * dense[left].y + alpha * dense[right].y});
    }
    return result;
  }

 private:
  std::vector<Point2> control_points_;
};

inline std::vector<double> solve_regularized_coordinate(
    const std::vector<Point2>& samples,
    int control_count,
    double smoothness,
    bool solve_x)
{
  const int unknown_count = control_count - 2;
  const int data_rows = static_cast<int>(samples.size());
  const int smoothing_rows = control_count - 2;
  const int row_count = data_rows + smoothing_rows;
  const int degree = control_count - 1;

  std::vector<double> matrix(
      static_cast<std::size_t>(row_count * unknown_count), 0.0);
  std::vector<double> rhs(
      static_cast<std::size_t>(std::max(row_count, unknown_count)), 0.0);

  const double first_value = solve_x ? samples.front().x : samples.front().y;
  const double last_value = solve_x ? samples.back().x : samples.back().y;

  for (int row = 0; row < data_rows; ++row) {
    const double t =
        static_cast<double>(row) / static_cast<double>(data_rows - 1);
    for (int control = 1; control < control_count - 1; ++control) {
      matrix[static_cast<std::size_t>(
          row * unknown_count + (control - 1))] =
          bernstein(control, degree, t);
    }
    const double value = solve_x ? samples[static_cast<std::size_t>(row)].x
                                 : samples[static_cast<std::size_t>(row)].y;
    rhs[static_cast<std::size_t>(row)] =
        value
        - first_value * bernstein(0, degree, t)
        - last_value * bernstein(degree, degree, t);
  }

  const double scale = std::sqrt(std::max(0.0, smoothness));
  for (int difference = 0; difference < control_count - 2; ++difference) {
    const int row = data_rows + difference;
    const int indices[3] = {difference, difference + 1, difference + 2};
    const double coefficients[3] = {1.0, -2.0, 1.0};
    double known_contribution = 0.0;

    for (int term = 0; term < 3; ++term) {
      const int control = indices[term];
      const double coefficient = scale * coefficients[term];
      if (control == 0) {
        known_contribution += coefficient * first_value;
      } else if (control == control_count - 1) {
        known_contribution += coefficient * last_value;
      } else {
        matrix[static_cast<std::size_t>(
            row * unknown_count + (control - 1))] += coefficient;
      }
    }
    rhs[static_cast<std::size_t>(row)] = -known_contribution;
  }

  const lapack_int status = LAPACKE_dgels(
      LAPACK_ROW_MAJOR,
      'N',
      static_cast<lapack_int>(row_count),
      static_cast<lapack_int>(unknown_count),
      1,
      matrix.data(),
      static_cast<lapack_int>(unknown_count),
      rhs.data(),
      1);
  if (status != 0) {
    throw std::runtime_error(
        "LAPACKE_dgels failed with status " + std::to_string(status));
  }

  rhs.resize(static_cast<std::size_t>(unknown_count));
  return rhs;
}

inline BezierCurve fit_regularized_bezier(
    const std::vector<Point2>& samples,
    int control_count,
    double smoothness)
{
  // The objective is the same regularized least-squares model used by
  // circleAirfoil.h::LSFittingSmooth():
  //
  //   min_P  sum_j |B(t_j; P) - X_j|^2
  //          + smoothness * sum_i |P_i - 2 P_{i+1} + P_{i+2}|^2,
  //
  // with the two end control points fixed at the surface endpoints.
  if (samples.size() < 3) {
    throw std::runtime_error("not enough samples for Bezier fitting");
  }
  if (control_count < 4
      || control_count >= static_cast<int>(samples.size()) + 2) {
    throw std::runtime_error("invalid number of Bezier control points");
  }

  const std::vector<double> internal_x =
      solve_regularized_coordinate(samples, control_count, smoothness, true);
  const std::vector<double> internal_y =
      solve_regularized_coordinate(samples, control_count, smoothness, false);

  std::vector<Point2> controls(static_cast<std::size_t>(control_count));
  controls.front() = samples.front();
  controls.back() = samples.back();
  for (int i = 0; i < control_count - 2; ++i) {
    controls[static_cast<std::size_t>(i + 1)] = {
        internal_x[static_cast<std::size_t>(i)],
        internal_y[static_cast<std::size_t>(i)]};
  }
  return BezierCurve(std::move(controls));
}

inline AirfoilData apply_gaussian_motion(
    const AirfoilData& original,
    double center,
    double width,
    double upper_shift,
    double lower_shift)
{
  if (!(width > 0.0)) {
    throw std::runtime_error("motion width must be positive");
  }

  AirfoilData moved = original;
  const auto displacement = [center, width](double x) {
    const double normalized = (x - center) / width;
    const double endpoint_window = std::max(0.0, 4.0 * x * (1.0 - x));
    return std::exp(-0.5 * normalized * normalized) * endpoint_window;
  };

  for (Point2& point : moved.upper) {
    point.y += upper_shift * displacement(point.x);
  }
  for (Point2& point : moved.lower) {
    point.y += lower_shift * displacement(point.x);
  }

  // Preserve the common leading edge exactly.
  const double leading_y =
      0.5 * (moved.upper.front().y + moved.lower.front().y);
  moved.upper.front().y = leading_y;
  moved.lower.front().y = leading_y;
  return moved;
}

struct FittedAirfoil {
  BezierCurve upper;
  BezierCurve lower;
  std::vector<Point2> boundary_points;
};

inline FittedAirfoil fit_and_sample(
    const AirfoilData& data,
    int control_count,
    double smoothness,
    int boundary_point_count)
{
  if (boundary_point_count < 12) {
    throw std::runtime_error("too few boundary points for EasyMesh");
  }

  BezierCurve upper =
      fit_regularized_bezier(data.upper, control_count, smoothness);
  BezierCurve lower =
      fit_regularized_bezier(data.lower, control_count, smoothness);

  const int upper_count = boundary_point_count / 2 + 1;
  const int lower_count = boundary_point_count - upper_count + 2;
  const std::vector<Point2> upper_points =
      upper.sample_uniform_arclength(upper_count);
  const std::vector<Point2> lower_points =
      lower.sample_uniform_arclength(lower_count);

  std::vector<Point2> boundary;
  boundary.reserve(static_cast<std::size_t>(boundary_point_count));

  // EasyMesh expects an internal boundary to be clockwise:
  // leading edge -> upper trailing edge -> lower trailing edge -> leading edge.
  boundary.insert(boundary.end(), upper_points.begin(), upper_points.end());
  // Both surfaces share the leading and trailing edges, so neither endpoint
  // is repeated when the lower surface is appended in reverse order.
  for (int i = lower_count - 2; i >= 1; --i) {
    boundary.push_back(lower_points[static_cast<std::size_t>(i)]);
  }

  if (static_cast<int>(boundary.size()) != boundary_point_count) {
    throw std::runtime_error("internal boundary point count mismatch");
  }
  return {std::move(upper), std::move(lower), std::move(boundary)};
}

inline void write_indexed_points(
    const std::filesystem::path& filename,
    const std::vector<Point2>& points)
{
  std::ofstream output(filename);
  if (!output) {
    throw std::runtime_error("cannot write points: " + filename.string());
  }
  output << std::setprecision(16) << std::fixed;
  for (std::size_t i = 0; i < points.size(); ++i) {
    output << i << ": " << points[i].x << " " << points[i].y << "\n";
  }
}

inline void write_xy_csv(
    const std::filesystem::path& filename,
    const std::vector<Point2>& points,
    const std::string& series)
{
  std::ofstream output(filename);
  if (!output) {
    throw std::runtime_error("cannot write CSV: " + filename.string());
  }
  output << "series,index,x,y\n";
  output << std::setprecision(16);
  for (std::size_t i = 0; i < points.size(); ++i) {
    output << series << "," << i << "," << points[i].x << "," << points[i].y
           << "\n";
  }
}

inline void write_fit_csv(
    const std::filesystem::path& filename,
    const AirfoilData& raw,
    const FittedAirfoil& fitted)
{
  std::ofstream output(filename);
  if (!output) {
    throw std::runtime_error("cannot write fit CSV: " + filename.string());
  }

  output << "kind,surface,index,x,y\n";
  output << std::setprecision(16);
  for (std::size_t i = 0; i < raw.upper.size(); ++i) {
    output << "raw,upper," << i << "," << raw.upper[i].x << ","
           << raw.upper[i].y << "\n";
  }
  for (std::size_t i = 0; i < raw.lower.size(); ++i) {
    output << "raw,lower," << i << "," << raw.lower[i].x << ","
           << raw.lower[i].y << "\n";
  }

  const std::vector<Point2> upper_curve =
      fitted.upper.sample_uniform_arclength(401);
  const std::vector<Point2> lower_curve =
      fitted.lower.sample_uniform_arclength(401);
  for (std::size_t i = 0; i < upper_curve.size(); ++i) {
    output << "curve,upper," << i << "," << upper_curve[i].x << ","
           << upper_curve[i].y << "\n";
  }
  for (std::size_t i = 0; i < lower_curve.size(); ++i) {
    output << "curve,lower," << i << "," << lower_curve[i].x << ","
           << lower_curve[i].y << "\n";
  }
  for (std::size_t i = 0;
       i < fitted.upper.control_points().size();
       ++i) {
    const Point2& point = fitted.upper.control_points()[i];
    output << "control,upper," << i << "," << point.x << "," << point.y
           << "\n";
  }
  for (std::size_t i = 0;
       i < fitted.lower.control_points().size();
       ++i) {
    const Point2& point = fitted.lower.control_points()[i];
    output << "control,lower," << i << "," << point.x << "," << point.y
           << "\n";
  }
  for (std::size_t i = 0; i < fitted.boundary_points.size(); ++i) {
    const Point2& point = fitted.boundary_points[i];
    output << "mesh-point,boundary," << i << "," << point.x << ","
           << point.y << "\n";
  }
}

inline double distance(const Point2& a, const Point2& b)
{
  return std::hypot(a.x - b.x, a.y - b.y);
}

inline void write_easymesh_domain(
    const std::filesystem::path& filename,
    const std::vector<Point2>& airfoil_points,
    double outer_center_x,
    double outer_center_y,
    double outer_radius,
    int outer_point_count,
    int outer_boundary_mark = 1,
    int airfoil_boundary_mark = 3)
{
  if (!(outer_radius > 1.0) || outer_point_count < 12) {
    throw std::runtime_error("invalid outer circle");
  }

  std::vector<Point2> outer(static_cast<std::size_t>(outer_point_count));
  for (int i = 0; i < outer_point_count; ++i) {
    const double theta =
        2.0 * std::acos(-1.0) * static_cast<double>(i)
        / static_cast<double>(outer_point_count);
    outer[static_cast<std::size_t>(i)] = {
        outer_center_x + outer_radius * std::cos(theta),
        outer_center_y + outer_radius * std::sin(theta)};
  }

  const int total_points =
      outer_point_count + static_cast<int>(airfoil_points.size());
  std::ofstream output(filename);
  if (!output) {
    throw std::runtime_error("cannot write EasyMesh input: " + filename.string());
  }
  output << std::setprecision(16) << std::fixed;
  output << total_points << "\n";

  int index = 0;
  for (int i = 0; i < outer_point_count; ++i) {
    const Point2& point = outer[static_cast<std::size_t>(i)];
    const Point2& next =
        outer[static_cast<std::size_t>((i + 1) % outer_point_count)];
    output << index++ << ": " << point.x << " " << point.y << " "
           << distance(point, next) << " " << outer_boundary_mark << "\n";
  }
  for (std::size_t i = 0; i < airfoil_points.size(); ++i) {
    const Point2& point = airfoil_points[i];
    const Point2& next = airfoil_points[(i + 1) % airfoil_points.size()];
    output << index++ << ": " << point.x << " " << point.y << " "
           << 10.0 * distance(point, next) << " "
           << airfoil_boundary_mark << "\n";
  }

  output << total_points << "\n";
  index = 0;
  for (int i = 0; i < outer_point_count; ++i) {
    output << index << ": " << index << " "
           << ((i + 1) % outer_point_count) << " "
           << outer_boundary_mark << "\n";
    ++index;
  }
  for (std::size_t i = 0; i < airfoil_points.size(); ++i) {
    const int local_next =
        static_cast<int>((i + 1) % airfoil_points.size());
    output << index << ": " << index << " "
           << outer_point_count + local_next << " "
           << airfoil_boundary_mark << "\n";
    ++index;
  }
}

}  // namespace airfoil_demo

#endif
