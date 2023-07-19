// @author Merve Asiler

#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/enum.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Quotient.h>
#include "FilteredPredicates.h"
#include "BasicGeometricElements.h"

typedef CGAL::Simple_cartesian<double>																	C;
typedef CGAL::Simple_cartesian<CGAL::Interval_nt_advanced>												InK;
typedef CGAL::Simple_cartesian<CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> >>					ExK;
typedef CGAL::Cartesian_converter<C, ExK>																C2ExK;
typedef CGAL::Cartesian_converter<C, InK>																C2InK;

typedef CGAL::Exact_predicates_inexact_constructions_kernel												IK;
typedef CGAL::Exact_kernel_selector<C>::Exact_kernel													Exact_kernel;
typedef CGAL::Exact_kernel_selector<C>::C2E																C2E;
typedef CGAL::Simple_cartesian<CGAL::Interval_nt_advanced>												Approximate_kernel;
typedef CGAL::Cartesian_converter<C, Approximate_kernel>												C2A;

int findClosestPlaneToPointin3D(vector<HalfSpace> halfSpaceSet, double* point);

int findClosestPlaneToPointin2D(vector<HalfSpace> halfSpaceSet, int base_plane_id, double* point);



// Filtered kernel to check whether the planes are parallel 
template <typename K>
struct ParallelPlanes {

	typedef typename K::Point_3					Point_3;
	typedef typename K::Plane_3					Plane_3;
	typedef typename bool						result_type;

	result_type operator()(	Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
							Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3) const
	{
		Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
		Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);

		if (CGAL::parallel(plane1, plane2))
			return true;
		return false;
	}

};

typedef CGAL::Filtered_predicate<ParallelPlanes<Exact_kernel>, ParallelPlanes<Approximate_kernel>, C2E, C2A>	PP;
//typedef CGAL::Filtered_predicate<ParallelPlanes<ExK>, ParallelPlanes<InK>, C2ExK, C2InK>	PP;

// Filtered kernel to check whether the planes are the same
template <typename K>
struct SamePlanes {

	typedef typename K::Point_3					Point_3;
	typedef typename K::Plane_3					Plane_3;
	typedef typename bool						result_type;

	result_type operator()(	Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
							Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3) const
	{
		Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
		Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);

		if (CGAL::parallel(plane1, plane2)) {
			if (CGAL::scalar_product(plane1.orthogonal_vector(), plane2.orthogonal_vector()) > 0)
				return true;
		}
		return false;

	}

};

typedef CGAL::Filtered_predicate<SamePlanes<Exact_kernel>, SamePlanes<Approximate_kernel>, C2E, C2A>	SP;
//typedef CGAL::Filtered_predicate<SamePlanes<ExK>, SamePlanes<InK>, C2ExK, C2InK>	SP;

// Filtered kernel to check whether the third plane gives an intersection with the first two 
template <typename K>
struct IntersectablePlane {

	typedef typename K::Point_3					Point_3;
	typedef typename K::Plane_3					Plane_3;
	typedef typename K::Line_3					Line_3;
	typedef typename CGAL::Sign					result_type;

	result_type operator()(	Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
							Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
							Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3) const
	{
		Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
		Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);

		Plane_3 candidate_plane(candidate_plane_p1, candidate_plane_p2, candidate_plane_p3);
		if (CGAL::parallel(candidate_plane, plane1) || CGAL::parallel(candidate_plane, plane2))
			return CGAL::NEGATIVE;	// not valid

		auto auto_line = CGAL::intersection(plane1, plane2);	// these two planes are assummed to be already intersectable
		Line_3* line = boost::get<Line_3>(&*auto_line);

		auto auto_candidate_point = CGAL::intersection(candidate_plane, *line);
		Point_3* candidate_point = boost::get<Point_3>(&*auto_candidate_point);

		if (candidate_point == NULL)
			return CGAL::NEGATIVE;	// not valid
		
		return CGAL::POSITIVE;

	}

};

typedef CGAL::Filtered_predicate<IntersectablePlane<Exact_kernel>, IntersectablePlane<Approximate_kernel>, C2E, C2A>	IP;
//typedef CGAL::Filtered_predicate<IntersectablePlane<ExK>, IntersectablePlane<InK>, C2ExK, C2InK>	IP;

// Filtered kernel to check whether the third plane gives an intersection with the first two and in front of the back one 
template <typename K>
struct IntersectableFrontPlane {

	typedef typename K::Point_3					Point_3;
	typedef typename K::Plane_3					Plane_3;
	typedef typename K::Line_3					Line_3;
	typedef typename CGAL::Sign					result_type;

	result_type operator()(	Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
							Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
							Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3,
							Point_3 back_plane_p1, Point_3 back_plane_p2, Point_3 back_plane_p3) const
	{
		Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
		Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);

		Plane_3 candidate_plane(candidate_plane_p1, candidate_plane_p2, candidate_plane_p3);
		if (CGAL::parallel(candidate_plane, plane1) || CGAL::parallel(candidate_plane, plane2))
			return CGAL::POSITIVE;	// not valid

		auto auto_line = CGAL::intersection(plane1, plane2);	// these two planes are assummed to be already intersectable
		Line_3* line = boost::get<Line_3>(&*auto_line);

		auto auto_candidate_point = CGAL::intersection(candidate_plane, *line);
		Point_3* candidate_point = boost::get<Point_3>(&*auto_candidate_point);

		if (candidate_point == NULL)
			return CGAL::POSITIVE;	// not valid

		CGAL::Sign sign = CGAL::orientation(back_plane_p1, back_plane_p2, back_plane_p3, *candidate_point);
		if (sign == CGAL::POSITIVE)
			return CGAL::POSITIVE;
		if (sign == CGAL::COPLANAR) {
			Plane_3 back_plane(back_plane_p1, back_plane_p2, back_plane_p3);
			if (CGAL::parallel(back_plane, candidate_plane) &&
				CGAL::scalar_product(back_plane.orthogonal_vector(), candidate_plane.orthogonal_vector()) > 0)
				return CGAL::POSITIVE;
		}
		return CGAL::NEGATIVE;
	}

};

typedef CGAL::Filtered_predicate<IntersectableFrontPlane<Exact_kernel>, IntersectableFrontPlane<Approximate_kernel>, C2E, C2A>	IFP;
//typedef CGAL::Filtered_predicate<IntersectableFrontPlane<ExK>, IntersectableFrontPlane<InK>, C2ExK, C2InK>	IFP;

// Filtered kernel to find the closest plane to a given point in 3D 
template <typename K>
struct ClosestPlaneToPoint_3D {
	typedef typename K::RT						RT;
	typedef typename K::Point_3					Point_3;
	typedef typename K::Plane_3					Plane_3;
	typedef typename CGAL::Sign					result_type;

	result_type operator()(	Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
							Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3, 
							Point_3 point) const
	{
		Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
		Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);

		RT distance1 = CGAL::squared_distance(plane1, point);
		RT distance2 = CGAL::squared_distance(plane2, point);
		return CGAL::compare(distance1, distance2);
	}

};

typedef CGAL::Filtered_predicate<ClosestPlaneToPoint_3D<Exact_kernel>, ClosestPlaneToPoint_3D<Approximate_kernel>, C2E, C2A>	CP_3D;
//typedef CGAL::Filtered_predicate<ClosestPlaneToPoint_3D<ExK>, ClosestPlaneToPoint_3D<InK>, C2ExK, C2InK>	CP_3D;

// Filtered kernel to find the closest plane to a given point on a plane in 2D 
template <typename K>
struct ClosestPlaneToPoint_2D {
	typedef typename K::RT						RT;
	typedef typename K::Point_3					Point_3;
	typedef typename K::Plane_3					Plane_3;
	typedef typename K::Line_3					Line_3;
	typedef typename CGAL::Sign					result_type;

	result_type operator()(	Point_3 base_plane_p1, Point_3 base_plane_p2, Point_3 base_plane_p3,
							Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
							Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
							Point_3 initialpoint) const
	{
		Plane_3 base_plane(base_plane_p1, base_plane_p2, base_plane_p3);
		Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
		Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);		// it is assummed that none of these planes are parallel

		Point_3 point = base_plane.projection(initialpoint);

		auto auto_line1 = CGAL::intersection(base_plane, plane1);
		Line_3* line1 = boost::get<Line_3>(&*auto_line1);

		auto auto_line2 = CGAL::intersection(base_plane, plane2);
		Line_3* line2 = boost::get<Line_3>(&*auto_line2);

		RT distance1 = CGAL::squared_distance(*line1, point);
		RT distance2 = CGAL::squared_distance(*line2, point);
		return CGAL::compare(distance1, distance2);
	}
};

typedef CGAL::Filtered_predicate<ClosestPlaneToPoint_2D<Exact_kernel>, ClosestPlaneToPoint_2D<Approximate_kernel>, C2E, C2A>	CP_2D;
//typedef CGAL::Filtered_predicate<ClosestPlaneToPoint_2D<ExK>, ClosestPlaneToPoint_2D<InK>, C2ExK, C2InK>	CP_2D;

// Filtered kernel to find the closest plane to a given point on a plane in 1D 
template <typename K>
struct ClosestPlaneToPoint_1D {

	typedef typename K::RT						RT;
	typedef typename K::Point_3					Point_3;
	typedef typename K::Plane_3					Plane_3;
	typedef typename K::Line_3					Line_3;
	typedef typename CGAL::Comparison_result	result_type;

	result_type operator()(	Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
							Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3, 
							Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3,
							Point_3 front_plane_p1, Point_3 front_plane_p2, Point_3 front_plane_p3,
							Point_3 back_plane_p1, Point_3 back_plane_p2, Point_3 back_plane_p3) const
	{
		Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
		Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);

		Plane_3 candidate_plane(candidate_plane_p1, candidate_plane_p2, candidate_plane_p3);
		if (CGAL::parallel(candidate_plane, plane1) || CGAL::parallel(candidate_plane, plane2))
			return CGAL::SMALLER;	// not valid

		auto auto_line = CGAL::intersection(plane1, plane2);
		Line_3* line = boost::get<Line_3>(&*auto_line);

		auto auto_candidate_point = CGAL::intersection(candidate_plane, *line);
		Point_3* candidate_point = boost::get<Point_3>(&*auto_candidate_point);

		if (candidate_point == NULL)
			return CGAL::SMALLER;

		CGAL::Sign sign = CGAL::orientation(back_plane_p1, back_plane_p2, back_plane_p3, *candidate_point);
		//cout << "SIGN: " << sign << endl;
		//Plane_3 back_plane2(back_plane_p1, back_plane_p2, back_plane_p3);
		//cout << "PLACEMENT: " << back_plane2.a() * candidate_point->x() + back_plane2.b() * candidate_point->y() + back_plane2.c() * candidate_point->z() + back_plane2.d() << endl;
		if (sign == CGAL::POSITIVE)
			return CGAL::SMALLER;
		if (sign == CGAL::COPLANAR) {
			Plane_3 back_plane(back_plane_p1, back_plane_p2, back_plane_p3);
			if (CGAL::parallel(back_plane, candidate_plane) &&
				CGAL::scalar_product(back_plane.orthogonal_vector(), candidate_plane.orthogonal_vector()) > 0)
				return CGAL::SMALLER;
		}

		Plane_3 front_plane(front_plane_p1, front_plane_p2, front_plane_p3);
		auto auto_front_point = CGAL::intersection(front_plane, *line);
		Point_3* front_point = boost::get<Point_3>(&*auto_front_point);

		Plane_3 back_plane(back_plane_p1, back_plane_p2, back_plane_p3);
		auto auto_back_point = CGAL::intersection(back_plane, *line);
		Point_3* back_point = boost::get<Point_3>(&*auto_back_point);
		
		RT distance1 = CGAL::squared_distance(*back_point, *front_point);
		RT distance2 = CGAL::squared_distance(*back_point, *candidate_point);
		return CGAL::compare(distance1, distance2);
	}

};

typedef CGAL::Filtered_predicate<ClosestPlaneToPoint_1D<Exact_kernel>, ClosestPlaneToPoint_1D<Approximate_kernel>, C2E, C2A>	CP_1D;
//typedef CGAL::Filtered_predicate<ClosestPlaneToPoint_1D<ExK>, ClosestPlaneToPoint_1D<InK>, C2ExK, C2InK>	CP_1D;

// Filtered kernel to find the closest plane to a given point on a plane in 1D 
template <typename K>
struct ClosestPlaneToPoint_UOD {

	typedef typename K::RT						RT;
	typedef typename K::Point_3					Point_3;
	typedef typename K::Plane_3					Plane_3;
	typedef typename K::Line_3					Line_3;
	typedef typename CGAL::Comparison_result	result_type;

	result_type operator()(	Point_3 plane1_p1, Point_3 plane1_p2, Point_3 plane1_p3,
							Point_3 plane2_p1, Point_3 plane2_p2, Point_3 plane2_p3,
							Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3,
							Point_3 front_plane_p1, Point_3 front_plane_p2, Point_3 front_plane_p3,
							Point_3 initialpoint) const
	{
		Plane_3 plane1(plane1_p1, plane1_p2, plane1_p3);
		Plane_3 plane2(plane2_p1, plane2_p2, plane2_p3);

		Plane_3 candidate_plane(candidate_plane_p1, candidate_plane_p2, candidate_plane_p3);
		if (CGAL::parallel(candidate_plane, plane1) || CGAL::parallel(candidate_plane, plane2))
			return CGAL::SMALLER;	// not valid

		auto auto_line = CGAL::intersection(plane1, plane2);
		Line_3* line = boost::get<Line_3>(&*auto_line);

		auto auto_candidate_point = CGAL::intersection(candidate_plane, *line);
		Point_3* candidate_point = boost::get<Point_3>(&*auto_candidate_point);

		if (candidate_point == NULL)
			return CGAL::SMALLER;

		Plane_3 front_plane(front_plane_p1, front_plane_p2, front_plane_p3);
		auto auto_front_point = CGAL::intersection(front_plane, *line);
		Point_3* front_point = boost::get<Point_3>(&*auto_front_point);

		Point_3 point = plane1.projection(initialpoint);
		point = line->projection(point);

		RT distance1 = CGAL::squared_distance(point, *front_point);
		RT distance2 = CGAL::squared_distance(point, *candidate_point);
		return CGAL::compare(distance1, distance2);
	}

};

typedef CGAL::Filtered_predicate<ClosestPlaneToPoint_UOD<Exact_kernel>, ClosestPlaneToPoint_UOD<Approximate_kernel>, C2E, C2A>	CP_UOD;
//typedef CGAL::Filtered_predicate<ClosestPlaneToPoint_UOD<ExK>, ClosestPlaneToPoint_UOD<InK>, C2ExK, C2InK>	CP_UOD;

// Filtered kernel to check whether the given plane is in front of a reference plane based on a reference point 
template <typename K>
struct FrontierPlane {

	typedef typename K::RT						RT;
	typedef typename K::Point_3					Point_3;
	typedef typename K::Plane_3					Plane_3;
	typedef typename K::Line_3					Line_3;
	typedef typename K::Vector_3				Vector_3;
	typedef typename CGAL::Sign					result_type;

	result_type operator()(	Point_3 base_plane_p1, Point_3 base_plane_p2, Point_3 base_plane_p3,
							Point_3 partner_plane_p1, Point_3 partner_plane_p2, Point_3 partner_plane_p3,
							Point_3 candidate_plane_p1, Point_3 candidate_plane_p2, Point_3 candidate_plane_p3,
							Point_3 reference_plane_p1, Point_3 reference_plane_p2, Point_3 reference_plane_p3) const
	{
		Plane_3 base_plane(base_plane_p1, base_plane_p2, base_plane_p3);
		Plane_3 partner_plane(partner_plane_p1, partner_plane_p2, partner_plane_p3);
		Plane_3 candidate_plane(candidate_plane_p1, candidate_plane_p2, candidate_plane_p3);
		Plane_3 reference_plane(reference_plane_p1, reference_plane_p2, reference_plane_p3);

		if (CGAL::parallel(candidate_plane, reference_plane))	// then they are the same planes
			return CGAL::POSITIVE;

		auto auto_line = CGAL::intersection(base_plane, candidate_plane);
		Line_3* line = boost::get<Line_3>(&*auto_line);
		auto auto_found_point = CGAL::intersection(*line, partner_plane);
		Point_3* found_point = boost::get<Point_3>(&*auto_found_point);

		Vector_3 direction = line->to_vector();
		Point_3 test_point = *found_point + direction;

		if (CGAL::orientation(partner_plane_p1, partner_plane_p2, partner_plane_p3, test_point) != CGAL::NEGATIVE) {
			direction = line->opposite().to_vector();
			test_point = *found_point + direction;
		}
		return CGAL::orientation(reference_plane_p1, reference_plane_p2, reference_plane_p3, test_point);
	}

};

typedef CGAL::Filtered_predicate<FrontierPlane<Exact_kernel>, FrontierPlane<Approximate_kernel>, C2E, C2A>	FP;
//typedef CGAL::Filtered_predicate<FrontierPlane<ExK>, FrontierPlane<InK>, C2ExK, C2InK>	FP;


