#include <cmath>
#include <cstdlib>
#include <iostream>

#include <AFEPack/AMGSolver.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Functional.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/Operator.h>
#include <AFEPack/TemplateElement.h>

namespace {

const double kPi = 4.0 * std::atan(1.0);

double exact_solution(const double* point) {
  return std::sin(kPi * point[0]) * std::sin(2.0 * kPi * point[1]);
}

double right_hand_side(const double* point) {
  return 5.0 * kPi * kPi * exact_solution(point);
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
              << " mesh_basename output_solution.dx\n";
    return EXIT_FAILURE;
  }

  EasyMesh mesh;
  mesh.readData(argv[1]);

  TemplateGeometry<2> triangle_template_geometry;
  triangle_template_geometry.readData("triangle.tmp_geo");
  CoordTransform<2, 2> triangle_coord_transform;
  triangle_coord_transform.readData("triangle.crd_trs");
  TemplateDOF<2> triangle_template_dof(triangle_template_geometry);
  triangle_template_dof.readData("triangle.1.tmp_dof");
  BasisFunctionAdmin<double, 2, 2> triangle_basis_function(
      triangle_template_dof);
  triangle_basis_function.readData("triangle.1.bas_fun");

  std::vector<TemplateElement<double, 2, 2> > template_element(1);
  template_element[0].reinit(triangle_template_geometry,
                             triangle_template_dof,
                             triangle_coord_transform,
                             triangle_basis_function);

  FEMSpace<double, 2> fem_space;
  fem_space.reinit(mesh, template_element);
  const int n_element = mesh.n_geometry(2);
  fem_space.element().resize(n_element);
  for (int i = 0; i < n_element; ++i) {
    fem_space.element(i).reinit(fem_space, i, 0);
  }
  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  StiffMatrix<2, double> stiffness(fem_space);
  stiffness.algebricAccuracy() = 4;
  stiffness.build();

  FEMFunction<double, 2> solution(fem_space);
  Vector<double> load;
  Operator::L2Discretize(&right_hand_side, fem_space, load, 4);

  BoundaryFunction<double, 2> boundary_1(
      BoundaryConditionInfo::DIRICHLET, 1, &exact_solution);
  BoundaryFunction<double, 2> boundary_2(
      BoundaryConditionInfo::DIRICHLET, 2, &exact_solution);
  BoundaryFunction<double, 2> boundary_3(
      BoundaryConditionInfo::DIRICHLET, 3, &exact_solution);
  BoundaryFunction<double, 2> boundary_4(
      BoundaryConditionInfo::DIRICHLET, 4, &exact_solution);
  BoundaryConditionAdmin<double, 2> boundary_admin(fem_space);
  boundary_admin.add(boundary_1);
  boundary_admin.add(boundary_2);
  boundary_admin.add(boundary_3);
  boundary_admin.add(boundary_4);
  boundary_admin.apply(stiffness, solution, load);

  AMGSolver solver(stiffness);
  solver.solve(solution, load, 1.0e-8, 200);
  solution.writeOpenDXData(argv[2]);

  const double error = Functional::L2Error(
      solution, FunctionFunction<double>(&exact_solution), 3);
  std::cout << "L2 error = " << error << '\n';
  return EXIT_SUCCESS;
}
