#include <deal.II/base/point.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <cmath>
#include <iostream>
#include <fstream>

#include <gtest/gtest.h>

using namespace dealii;

//! Generate a hypercube, and output it as an svg file.
void first_grid(Triangulation<2> &triangulation) {
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(4);

  std::ofstream out("grid-1.svg");
  GridOut grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to grid-1.svg" << std::endl;
}

//! Generate a locally refined hyper_shell, and output it as an svg file.
void second_grid(Triangulation<2> &triangulation) {
  const Point<2> center(1, 0);
  const double inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(triangulation, center, inner_radius, outer_radius,
                             10);

  // triangulation.reset_manifold(0);

  for (unsigned int step = 0; step < 5; ++step) {
    for (auto &cell : triangulation.active_cell_iterators()) {
      for (const auto v : cell->vertex_indices()) {
        const double distance_from_center = center.distance(cell->vertex(v));

        if (std::fabs(distance_from_center - inner_radius) <=
            1e-6 * inner_radius) {
          cell->set_refine_flag();
          break;
        }
      }
    }

    triangulation.execute_coarsening_and_refinement();
  }

  std::ofstream out("grid-2.svg");
  GridOut grid_out;
  grid_out.write_svg(triangulation, out);

  std::cout << "Grid written to grid-2.svg" << std::endl;
}

//! Create an L-shaped domain with one global refinement, and write it on
// `third_grid.vtk`.  Refine the L-shaped mesh adaptively around the re-entrant
// corner three times (after the global refinement you already did), but with a
// twist: refine all cells with the distance between the center of the cell and
// re-entrant corner is smaller than 1/3.
void third_grid(Triangulation<2> &tria) {
  // Insert code here
}

//! Returns a tuple with number of levels, number of cells, number of active
// cells. Test this with all of  your meshes.
std::tuple<unsigned int, unsigned int, unsigned int>
get_info(const Triangulation<2> &) {
  // Insert code here
  return std::make_tuple(0, 0, 0);
}

// int main() {
//   {
//     Triangulation<2> triangulation;
//     first_grid(triangulation);
//   }
//   {
//     Triangulation<2> triangulation;
//     second_grid(triangulation);
//   }
// }

TEST(Step1, Mark1) {
  Triangulation<2> tria;
  first_grid(tria);
  ASSERT_TRUE(std::ifstream("grid-1.svg"));
}

TEST(Step1, Mark2) {
  Triangulation<2> tria;
  second_grid(tria);
  ASSERT_TRUE(std::ifstream("grid-2.svg"));
}

TEST(Step1, Mark3) {
  Triangulation<2> tria;
  third_grid(tria);
  ASSERT_TRUE(std::ifstream("grid-3.vtk"));
}

TEST(Step1, Mark4) {
  Triangulation<2> tria;
  first_grid(tria);
  auto [levels, cells, active_cells] = get_info(tria);
  EXPECT_EQ(levels, 5u);
  EXPECT_EQ(cells, 341u);
  EXPECT_EQ(active_cells, 256u);
}

TEST(Step1, Mark5) {
  Triangulation<2> tria;
  second_grid(tria);
  const auto [levels, cells, active_cells] = get_info(tria);
  EXPECT_EQ(levels, 6u);
  EXPECT_EQ(cells, 1250u);
  EXPECT_EQ(active_cells, 940u);
}

TEST(Step1, Mark6) {
  Triangulation<2> tria;
  third_grid(tria);
  auto [levels, cells, active_cells] = get_info(tria);
  EXPECT_EQ(levels, 5u);
  EXPECT_EQ(cells, 351u);
  EXPECT_EQ(active_cells, 264u);
}
