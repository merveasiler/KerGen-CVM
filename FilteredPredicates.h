// @author Merve Asiler

#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/enum.h>

typedef CGAL::Simple_cartesian<double>							C;
typedef CGAL::Exact_predicates_inexact_constructions_kernel		IK;
typedef CGAL::Exact_predicates_exact_constructions_kernel		EK;
typedef CGAL::Cartesian_converter<C, EK>						C2EK;
typedef CGAL::Cartesian_converter<C, IK>						C2IK;
typedef CGAL::Cartesian_converter<EK, C>						E2C;
typedef CGAL::Cartesian_converter<EK, IK>						E2I;

/* Filtered kernel to find the closest plane to a given point in 3D */
template <typename K>
struct closest_plane_to_point_in_3D {
	typedef typename K::RT						RT;
	typedef typename K::Point_3					Point_3;
	typedef typename K::Plane_3					Plane_3;
	typedef typename CGAL::Sign					result_type;

	result_type operator()(const Point_3& point,
		const Point_3& closest_plane_point1,
		const Point_3& closest_plane_point2,
		const Point_3& closest_plane_point3,
		const Point_3& new_plane_point1,
		const Point_3& new_plane_point2,
		const Point_3& new_plane_point3) const
	{
		Plane_3 closest_plane(closest_plane_point1, closest_plane_point2, closest_plane_point3);
		Plane_3 new_plane(new_plane_point1, new_plane_point2, new_plane_point3);
		Point_3 closest_point = closest_plane.projection(point);
		Point_3 new_point = new_plane.projection(point);
		return CGAL::compare_distance_to_point(point, closest_point, new_point);
	}
};

typedef CGAL::Filtered_predicate<closest_plane_to_point_in_3D<EK>, closest_plane_to_point_in_3D<IK>, C2EK, C2IK>	CPP_3D;



/* Filtered kernel to find the closest plane to a given point on a plane in 2D */
template <typename K>
struct closest_plane_to_point_in_2D {
	typedef typename K::RT						RT;
	typedef typename K::Point_3					Point_3;
	typedef typename K::Plane_3					Plane_3;
	typedef typename K::Line_3					Line_3;
	typedef typename K::Compare_distance_3		cmp;
	typedef typename CGAL::Sign					result_type;

	result_type operator()(const Point_3& point,
		const Point_3& base_plane_point1,
		const Point_3& base_plane_point2,
		const Point_3& base_plane_point3,
		const Point_3& closest_plane_point1,
		const Point_3& closest_plane_point2,
		const Point_3& closest_plane_point3,
		const Point_3& new_plane_point1,
		const Point_3& new_plane_point2,
		const Point_3& new_plane_point3) const
	{
		Plane_3 base_plane(base_plane_point1, base_plane_point2, base_plane_point3);
		Plane_3 closest_plane(closest_plane_point1, closest_plane_point2, closest_plane_point3);
		Plane_3 new_plane(new_plane_point1, new_plane_point2, new_plane_point3);

		auto line1 = CGAL::intersection(base_plane, closest_plane);
		Line_3* closest_line = boost::get<Line_3>(&*line1);

		auto line2 = CGAL::intersection(base_plane, new_plane);
		Line_3* new_line = boost::get<Line_3>(&*line2);

		Point_3 closest_point = closest_line->projection(point);
		Point_3 new_point = new_line->projection(point);
		return CGAL::compare_distance_to_point(point, closest_point, new_point);

	}
};

typedef CGAL::Filtered_predicate<closest_plane_to_point_in_2D<EK>, closest_plane_to_point_in_2D<IK>, C2EK, C2IK>		CPP_2D;

/* Filtered kernel to find the closest plane to a given point on a line in 2D */
template <typename K>
struct closest_plane_to_point_on_line {
	typedef typename K::RT						RT;
	typedef typename K::Point_3					Point_3;
	typedef typename K::Plane_3					Plane_3;
	typedef typename K::Vector_3				Vector_3;
	typedef typename K::Ray_3					Ray_3;
	typedef typename CGAL::Sign					result_type;

	result_type operator()( const Point_3& point,
							const Point_3& ray_point,
							const Point_3& ray_plane1_point1,
							const Point_3& ray_plane1_point2,
							const Point_3& ray_plane1_point3,
							const Point_3& ray_plane2_point1,
							const Point_3& ray_plane2_point2,
							const Point_3& ray_plane2_point3,
							const Point_3& closest_plane_point1,
							const Point_3& closest_plane_point2,
							const Point_3& closest_plane_point3,
							const Point_3& new_plane_point1,
							const Point_3& new_plane_point2,
							const Point_3& new_plane_point3) const
	{

		Plane_3 ray_plane1(ray_plane1_point1, ray_plane1_point2, ray_plane1_point3);
		Plane_3 ray_plane2(ray_plane2_point1, ray_plane2_point2, ray_plane2_point3);
		Vector_3 ray_vector = CGAL::cross_product(ray_plane1.orthogonal_vector(), ray_plane2.orthogonal_vector());
		Ray_3 ray(ray_point, ray_vector);
		Plane_3 closest_plane(closest_plane_point1, closest_plane_point2, closest_plane_point3);
		Plane_3 new_plane(new_plane_point1, new_plane_point2, new_plane_point3);

		auto point1 = CGAL::intersection(closest_plane, ray);
		Point_3* closest_point = boost::get<Point_3>(&*point1);

		auto point2 = CGAL::intersection(new_plane, ray);
		if (!point2)
			return CGAL::SMALLER;
		Point_3* new_point = boost::get<Point_3>(&*point2);

		return CGAL::compare_distance_to_point(point, *closest_point, *new_point);
	}
};

typedef CGAL::Filtered_predicate<closest_plane_to_point_on_line<EK>, closest_plane_to_point_on_line<IK>, C2EK, C2IK>		CPP_OnLine;

/* Filtered kernel to find the right index for ordering the planes */
/*
template <typename K>	
struct right_index_for_the_plane {
	typedef typename K::RT						RT;
	typedef typename K::Point_3					Point_3;
	typedef typename K::Plane_3					Plane_3;
	typedef typename K::Vector_3				Vector_3;
	typedef typename K::Ray_3					Ray_3;
	typedef typename CGAL::Sign					result_type;

	result_type operator()(	const Point_3& point,
							const Point_3& ray_point,
							vector<const Point_3>& plane_points,
							vector<vector<int>>& plane_point_ids,
							int base_plane_id, 
							int partner_plane_id, 
							vector<int>& next_partner_ids) const
	{

		vector<int> ordered_next_partner_ids;
		vector<EK::Vector_3> array_of_direction_between_partners;
		vector<EK::Vector_3> array_of_next_edge_direction;

		vector<Plane_3> next_partner_planes;
		for (int i = 0; i < next_partner_ids.size(); i++)
			planes.push_back(Plane_3(	plane_points[next_partner_ids[i][0]], 
										plane_points[next_partner_ids[i][1]], 
										plane_points[next_partner_ids[i][2]])	);
		Plane_3 base_plane = Plane_3(	plane_points[base_plane_id][0],
										plane_points[base_plane_id][1],
										plane_points[base_plane_id][2]);
		Plane_3 partner_plane = Plane_3(plane_points[partner_plane_id][0],
										plane_points[partner_plane_id][1],
										plane_points[partner_plane_id][2]);

		Vector_3 incoming_ray_direction = CGAL::cross_product(	base_plane.orthogonal_vector(), partner_plane.orthogonal_vector());

		// to order...
		for (int i = 0; i < next_partner_ids.size(); i++) {

			Vector_3 direction_between_partners = CGAL::cross_product(partner_plane.orthogonal_vector(), next_partner_planes[i].orthogonal_vector());

			Vector_3 next_edge_direction = CGAL::cross_product(base_plane.orthogonal_vector(), next_partner_planes[i].orthogonal_vector());

			Ray_3 ray1(point, direction_between_partners);
		

			// first check book page case
			bool book_page_detected = false;
			bool no_change_then_continue = false;
			for (int j = 0; j < ordered_next_partner_ids.size(); j++) {
				Ray_3 ray2(point, array_of_direction_between_partners[j]);
				if (CGAL::parallel(ray1, ray2)) {	// we should eliminate one
					if (CGAL::compare_dihedral_angle(array_of_direction_between_partners[j], incoming_ray_direction, array_of_next_edge_direction[j],
						direction_between_partners, incoming_ray_direction, next_edge_direction) == CGAL::LARGER) {
						// eliminate the previous j^th one
						int k = j;
						for (; k < ordered_next_partner_ids.size() - 1; k++) {
							array_of_direction_between_partners[k] = array_of_direction_between_partners[k+1];
							array_of_next_edge_direction[k] = array_of_next_edge_direction[k+1];
							ordered_next_partner_ids[k] = ordered_next_partner_ids[k+1];
						}
						array_of_direction_between_partners[k] = direction_between_partners;
						array_of_next_edge_direction[k] = next_edge_direction;
						ordered_next_partner_ids[k] = next_partner_ids[i];
					}
					else	// do not take the current one (i^th) into account
						no_change_then_continue = true;
					book_page_detected = true;
					break;
				}
			}

			if (no_change_then_continue)
				continue;
			if (!book_page_detected) {
				array_of_direction_between_partners.push_back(direction_between_partners);
				array_of_next_edge_direction.push_back(next_edge_direction);
				ordered_next_partner_ids.push_back(next_partner_ids[i]);
			}

			// then order
			int j = ordered_next_partner_ids.size() - 2;
			for ( ; j >= 0; j--) {
				// the one whose dihedral angle is smaller should be processed later
				if (CGAL::compare_dihedral_angle(base_plane.orthogonal_vector(), incoming_ray_direction, array_of_next_edge_direction[j],
					base_plane.orthogonal_vector(), incoming_ray_direction, next_edge_direction) == CGAL::SMALLER) {
					array_of_direction_between_partners[j + 1] = array_of_direction_between_partners[j];
					array_of_next_edge_direction[j + 1] = array_of_next_edge_direction[j];
					ordered_next_partner_ids[j + 1] = ordered_next_partner_ids[j];
					continue;
				}
				else
					break;
			}
			array_of_direction_between_partners[j + 1] = direction_between_partners;
			array_of_next_edge_direction[j + 1] = next_edge_direction;
			ordered_next_partner_ids[j + 1] = next_partner_ids[i];
			break;

		}


	}
};

typedef CGAL::Filtered_predicate<right_index_for_the_plane<EK>, right_index_for_the_plane<IK>, C2EK, C2IK>		RIP;
*/