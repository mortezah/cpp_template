#ifndef CGAL_TYPEDEFS_H
#define CGAL_TYPEDEFS_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <math.h>

#include <expected>
#include <format>
#include <optional>
#include <random>

// kernels and traits typedefs
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Surface_mesh = CGAL::Surface_mesh<Point_3>;
/* using Filtered_graph = CGAL::Face_filtered_graph<Surface_mesh>; */

// boost graph-traits related typedefs
using Graph_traits = boost::graph_traits<Surface_mesh>;
using vertex_iterator = Graph_traits::vertex_iterator;
using face_iterator = Graph_traits::face_iterator;
using halfedge_iterator = Graph_traits::halfedge_iterator;

// entity descriptor typedefs
using vertex_descriptor = Graph_traits::vertex_descriptor;
using face_descriptor = Graph_traits::face_descriptor;
using edge_descriptor = Graph_traits::edge_descriptor;
using halfedge_descriptor = Graph_traits::halfedge_descriptor;

/* // entity size types */
using faces_stype = Graph_traits::faces_size_type;

#endif  // CGAL_TYPEDEFS_H
