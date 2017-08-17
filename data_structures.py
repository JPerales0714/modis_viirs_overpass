from math import radians, sin, cos, asin, sqrt, atan2, degrees
import numpy


class NadirPoint(object):

    def __init__(self, lat, long, scan_time, scan_number, spatial_search_range, temporal_search_range, nadir_along_frame_index):
        self.coordinate = (lat,long)
        self.scan_time = scan_time
        self.scan_number = scan_number
        self.max_coordinate_difference = spatial_search_range
        self.max_time_difference = temporal_search_range
        self.nadir_swath_frame = nadir_along_frame_index

    def __str__(self):
        return "scan: " + str(self.scan_number) + " location: " + str(self.coordinate) + " scan time: " + str(self.scan_time)

    def get_coordinates(self):
        return self.coordinate

    def get_nadir_pos(self):
        return self.nadir_swath_frame

    def get_time(self):
        return self.scan_time

    def within_time_range(self, other_time):
        time_difference = self.scan_time - other_time
        if time_difference <= self.max_time_difference:
            return True
        else:
            return False

    def within_geospatial_range(self, other_coordinate):
        latitude_difference = self.coordinate[0] - other_coordinate[0]
        longitude_difference = self.coordinate[1] - other_coordinate[1]
        if abs(latitude_difference) <= self.max_coordinate_difference:
            if abs(longitude_difference) <= (self.max_coordinate_difference + self.max_coordinate_difference*(1-cos(self.coordinate[0]))):
                return True
            else:
                return False
        else:
            return False

class TwoPointComparison(object):
    # TODO: Improve Constructor Clarity

    def __init__(self, viirs_scan_number, modis_scan_number, viirs_swath_position, modis_swath_position, modis_loc=(0, 0), viirs_loc=(0, 0)):
        self.viirs_scan_number = viirs_scan_number
        self.modis_scan_number = modis_scan_number
        self.viirs_swath_pos = viirs_swath_position
        self.modis_swath_pos = modis_swath_position
        self.viirs_location = viirs_loc
        self.modis_location = modis_loc
        self.viirs_value = 0.0
        self.modis_value = 0.0
        self.difference_ratio = 0.0
        self.scan_angle = 0.0

    def set_comparison_values_viirs_offnad(self, v_rad, m_rad, swath):
        self.viirs_value = v_rad
        self.modis_value = m_rad
        self.difference_ratio = v_rad / m_rad
        self.scan_angle = (2.*swath +0.5-33.5)*0.017785 -56.063

    def set_comparison_values_modis_offnad(self, v_rad, m_rad, swath):
        self.viirs_value = v_rad
        self.modis_value = m_rad
        self.difference_ratio = m_rad / v_rad
        self.scan_angle = 2.0 * ((10.5+swath/1353 * 55.0) - 38.0)

    def return_info(self):
        tp = (self.viirs_scan_number, self.modis_scan_number, self.viirs_swath_pos, self.modis_swath_pos, self.viirs_value, self.modis_value, self.difference_ratio, self.scan_angle)
        return tp

    def get_ratio(self):
        return self.difference_ratio

    def get_angle(self):
        return self.scan_angle

    def __str__(self):
        return "MODIS SCAN: " + str(self.modis_scan_number) + " MODIS RAD: " + str(self.modis_value) + \
               " VIIRS SCAN: " + str(self.viirs_scan_number) + " VIIRS RAD: " + str(self.viirs_value)

    def get_viirs_scan(self):
        return self.viirs_scan_number

    def get_modis_scan(self):
        return self.modis_scan_number

    def get_viirs_swath_pos(self):
        return self.viirs_swath_pos

    def get_modis_swath_pos(self):
        return self.modis_swath_pos


class GeospatialScanBox(object):

    def __init__(self, tl_coordinate, bl_coordinate, tr_coordinate, br_coordinate, start_scan, end_scan):
        self.top_left = tl_coordinate
        self.top_right = tr_coordinate
        self.bottom_right = br_coordinate
        self.bottom_left = bl_coordinate
        self.start_scan = start_scan
        self.end_scan = end_scan
        self.crosses_anti = self.crosses_antemeridian()

    def crosses_antemeridian(self):
        if abs(self.top_left[1]) > 80.0 and abs(self.top_right[1]) > 80.0:
            return not ((self.top_left[1] >= 0.0 and self.top_right[1] >= 0.0) or (self.top_left[1] < 0.0 and self.top_right[1] < 0.0))
        else:
            return False

    def get_starting_scan(self):
        return self.start_scan

    def get_ending_scan(self):
        return self.end_scan

    def get_brng_d(self, lat1, lon1, lat2, lon2):
        if self.crosses_anti:
            if (-180 <= lon1 < 0):
                lon1 += 360
            if (-180 <= lon2 < 0):
                lon2 += 360
        lat1 = radians(lat1)
        lon1 = radians(lon1)
        lat2 = radians(lat2)
        lon2 = radians(lon2)

        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
        c = 2 * asin(sqrt(a))
        d = c * 6371
        y = sin(dlon) * cos(lat2)
        x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon)
        brng = atan2(y, x)
        return d, brng

    def get_new_point(self, lat1, lon1, d, brng):
        if self.crosses_anti and (-180 <= lon1 < 0):
            lon1 += 360
        R = 6371
        lat1 = radians(lat1)
        lon1 = radians(lon1)
        lat2 = asin(sin(lat1) * cos(d / R) + cos(lat1) * sin(d / R) * cos(brng))
        lon2 = lon1 + atan2(sin(brng) * sin(d / R) * cos(lat1), cos(d / R) - sin(lat1) * sin(lat2))

        lat2 = degrees(lat2)
        lon2 = degrees(lon2)
        if self.crosses_anti and lon2 > 180:
            lon2 -= 360
        return lat2, lon2

    def get_edge_points(self, point1, point2):
        distance, bearing = self.get_brng_d(point1[0], point1[1], point2[0], point2[1])
        increment = distance / 30
        point_list = []
        for i in numpy.arange(0, distance, increment):
            point_list.append(tuple(self.get_new_point(point1[0], point1[1], i, bearing)))
        return point_list

    def get_edge_values(self):
        edge1 = self.get_edge_points(self.top_left, self.top_right)
        edge2 = self.get_edge_points(self.top_right, self.bottom_right)
        edge3 = self.get_edge_points(self.bottom_left, self.bottom_right)
        edge4 = self.get_edge_points(self.top_left, self.bottom_left)
        return edge1 + edge2 + edge3 + edge4

    def encapsulates(self, coordinate):
        edge_points = self.get_edge_values()
        smallest_lat = 1000
        smallest_lon = 1000
        biggest_lat = -1000
        biggest_lon = -1000
        for point in edge_points:
            if point[0] > biggest_lat:
                biggest_lat = point[0]
            if point[0] < smallest_lat:
                smallest_lat = point[0]
            if point[1] > biggest_lon:
                biggest_lon = point[1]
            if point[1] < smallest_lon:
                smallest_lon = point[1]
        if (smallest_lon <= coordinate[1] <= biggest_lon) and (smallest_lat <= coordinate[0] <= biggest_lat):
            return True
        else:
            return False

    def get_scan_list(self):
        list = []
        # end_scan + 1 because the right bound is not inclusive
        for i in range(self.start_scan, self.end_scan+1):
            list.append(i)
        return list

    def __str__(self):
        return str(self.top_left) + " " + str(self.bottom_right)
