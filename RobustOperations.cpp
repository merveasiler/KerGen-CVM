#include "RobustOperations.h"

IK::Plane_3 defineInexactPlane(HalfSpace halfspace) {

	IK::Plane_3 plane(	IK::Point_3(halfspace.point1[0], halfspace.point1[1], halfspace.point1[2]),
						IK::Point_3(halfspace.point2[0], halfspace.point2[1], halfspace.point2[2]),
						IK::Point_3(halfspace.point3[0], halfspace.point3[1], halfspace.point3[2]));
	return plane;

}

EK::Plane_3 defineExactPlane(HalfSpace halfspace) {

	EK::Plane_3 plane(	EK::Point_3(halfspace.point1[0], halfspace.point1[1], halfspace.point1[2]),
						EK::Point_3(halfspace.point2[0], halfspace.point2[1], halfspace.point2[2]),
						EK::Point_3(halfspace.point3[0], halfspace.point3[1], halfspace.point3[2]));
	return plane;

}

IK::Point_3 defineInexactPoint(double* point) {

	IK::Point_3 ik_point(point[0], point[1], point[2]);
	return ik_point;

}

EK::Point_3 defineExactPoint(double* point) {

	EK::Point_3 ek_point(point[0], point[1], point[2]);
	return ek_point;

}

int findClosestPlaneToPointin3D(vector<HalfSpace> halfSpaceSet, double* point) {

	/*  exceptions related to inexact computations
	* 
		catch (CGAL::Precondition_exception& e) {
			// Handle precondition exception
		}
		catch (CGAL::Assertion_exception& e) {
			// Handle assertion exception
		}
		catch (CGAL::Inexact_conversion_exception& e) {
			// Handle inexact conversion exception
		}
		catch (CGAL::Null_result_exception& e) {
			// Handle null result exception
		}
	*/

	IK::Point_3 ik_point = defineInexactPoint(point);
	EK::Point_3 ek_point = defineExactPoint(point);

	// initialize the closest plane
	int closest_plane_id = 0;
	CGAL::Lazy_exact_nt<boost::multiprecision::mpq_rational> ik_closest_distance;
	try {
		IK::Plane_3 ik_closest_plane = defineInexactPlane(halfSpaceSet[closest_plane_id]);
		ik_closest_distance = CGAL::squared_distance(ik_closest_plane, ik_point);
	} 
	catch (...) {
		EK::Plane_3 ek_closest_plane = defineExactPlane(halfSpaceSet[closest_plane_id]);
		ik_closest_distance = CGAL::squared_distance(ek_closest_plane, ek_point);
	}
	
	// check for the other planes if they are closer
	for (int i = 1; i < halfSpaceSet.size(); i++) {

		int current_plane_id = i;

		try {

			IK::Plane_3 ik_current_plane = defineInexactPlane(halfSpaceSet[current_plane_id]);
			CGAL::Lazy_exact_nt<boost::multiprecision::mpq_rational> ik_current_distance = CGAL::squared_distance(ik_current_plane, ik_point);
			
			if (ik_current_distance < ik_closest_distance) {
				closest_plane_id = current_plane_id;
				ik_closest_distance = ik_current_distance;
			}

		}
		catch (...) {
			EK::Plane_3 ek_current_plane = defineExactPlane(halfSpaceSet[current_plane_id]);
			CGAL::Lazy_exact_nt<boost::multiprecision::mpq_rational> ek_current_distance = CGAL::squared_distance(ek_current_plane, ek_point);

			EK::Plane_3 ek_closest_plane = defineExactPlane(halfSpaceSet[closest_plane_id]);
			CGAL::Lazy_exact_nt<boost::multiprecision::mpq_rational> ek_closest_distance = CGAL::squared_distance(ek_closest_plane, ek_point);

			if (ek_current_distance < ek_closest_distance) {
				closest_plane_id = current_plane_id;
				ik_closest_distance = ek_current_distance;
			}

		}

	}

	return closest_plane_id;
}



int findClosestPlaneToPointin2D(vector<HalfSpace> halfSpaceSet, int base_plane_id, double* point) {

	EK::Plane_3 ek_base_plane = defineExactPlane(halfSpaceSet[base_plane_id]);
	IK::Plane_3 ik_base_plane = defineInexactPlane(halfSpaceSet[base_plane_id]);
	
	EK::Point_3 ek_point = ek_base_plane.projection(defineExactPoint(point));
	IK::Point_3 ik_point = ik_base_plane.projection(defineInexactPoint(point));

	// initialize the closest plane
	CGAL::Lazy_exact_nt<boost::multiprecision::mpq_rational> ik_closest_distance;
	int closest_plane_id = -1;

	for (int i = 0; i < halfSpaceSet.size(); i++) {

		if (i == base_plane_id)
			continue;

		int current_plane_id = i;
		CGAL::Lazy_exact_nt<boost::multiprecision::mpq_rational> ik_current_distance;
		try {
			IK::Plane_3 ik_current_plane = defineInexactPlane(halfSpaceSet[current_plane_id]);
			if (CGAL::parallel(ik_base_plane, ik_current_plane))
				throw - 1;
			auto auto_common_line = CGAL::intersection(ik_base_plane, ik_current_plane);
			IK::Line_3* common_line = boost::get<IK::Line_3>(&*auto_common_line);
			ik_current_distance = CGAL::squared_distance(*common_line, ik_point);
		}
		catch(...) {
			EK::Plane_3 ek_current_plane = defineExactPlane(halfSpaceSet[current_plane_id]);
			if (CGAL::parallel(ek_base_plane, ek_current_plane))
				continue;
			auto auto_common_line = CGAL::intersection(ek_base_plane, ek_current_plane);
			EK::Line_3* common_line = boost::get<EK::Line_3>(&*auto_common_line);
			ik_current_distance = CGAL::squared_distance(*common_line, ek_point);
		}

		closest_plane_id = current_plane_id;
		ik_closest_distance = ik_current_distance;
		break;
	}

	// check for the other planes if they are closer
	for (int i = closest_plane_id + 1; i < halfSpaceSet.size(); i++) {

		if (i == base_plane_id)
			continue;

		int current_plane_id = i;
		try {
			IK::Plane_3 ik_current_plane = defineInexactPlane(halfSpaceSet[current_plane_id]);
			if (CGAL::parallel(ik_base_plane, ik_current_plane))
				throw - 1;
			auto auto_common_line = CGAL::intersection(ik_base_plane, ik_current_plane);
			IK::Line_3* common_line = boost::get<IK::Line_3>(&*auto_common_line);
			CGAL::Lazy_exact_nt<boost::multiprecision::mpq_rational> ik_current_distance = CGAL::squared_distance(*common_line, ik_point);

			if (ik_current_distance < ik_closest_distance) {
				closest_plane_id = current_plane_id;
				ik_closest_distance = ik_current_distance;
			}
		}

		catch (int e) {
			cout << "!!! EXCEPTION !!! " << e << endl;
			EK::Plane_3 ek_current_plane = defineExactPlane(halfSpaceSet[current_plane_id]);
			if (CGAL::parallel(ek_base_plane, ek_current_plane))
				continue;
			auto auto_common_line = CGAL::intersection(ek_base_plane, ek_current_plane);
			EK::Line_3* common_line = boost::get<EK::Line_3>(&*auto_common_line);
			CGAL::Lazy_exact_nt<boost::multiprecision::mpq_rational> ek_current_distance = CGAL::squared_distance(*common_line, ek_point);

			EK::Plane_3 ek_closest_plane = defineExactPlane(halfSpaceSet[closest_plane_id]);
			auto auto_closest_common_line = CGAL::intersection(ek_base_plane, ek_closest_plane);
			EK::Line_3* closest_common_line = boost::get<EK::Line_3>(&*auto_closest_common_line);
			CGAL::Lazy_exact_nt<boost::multiprecision::mpq_rational> ek_closest_distance = CGAL::squared_distance(*closest_common_line, ek_point);

			if (ek_current_distance < ek_closest_distance) {
				closest_plane_id = current_plane_id;
				ik_closest_distance = ek_current_distance;
			}
		}

	}

	return closest_plane_id;
}