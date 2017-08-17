from pyhdf.SD import *
from pyhdf.HDF import *
import h5py
import re
import time
import numpy
import datetime
from data_sets import AquaSDSDataSet, SuomiDataSet, AquaVDataSet
from data_structures import NadirPoint, TwoPointComparison, GeospatialScanBox


class HDF4File(object):

    def __init__(self, file_name):
        self.hdf_file = HDF(file_name, HC.READ)
        self.sd_file_interface = SD(file_name, SDC.READ)
        self.v_file_interface = self.hdf_file.vstart()
        self.attributes = self.sd_file_interface.attributes()
        self.boxes = self.generate_lat_lon_boxes()
        self.name = file_name.split("/")[-1]

    def get_attributes(self):
        return self.attributes

    def get_times_list(self):
        start_time = self.get_start_time()
        scan_number = self.get_number_of_scans()
        times = []
        for i in range(scan_number):
            times.append(start_time + datetime.timedelta(seconds=(1.4771 * i)))
        return times

    def get_start_time(self):
        core_metadata = self.attributes['CoreMetadata.0']
        regex = re.compile('[M][Y][D][.A-Z0-9]+[.][h][d][f]')
        file_descriptor = re.search(regex, core_metadata).group(0).split(".")
        # Fetching date information
        year = int(file_descriptor[1][1:5])
        day_of_year = int(file_descriptor[1][5:])
        date = datetime.datetime(year, 1, 1) + datetime.timedelta(days=day_of_year-1)
        # Fetching time information
        start_hour = int(file_descriptor[2][0:2])
        start_minute = int(file_descriptor[2][-2:])
        date_and_time = date + datetime.timedelta(hours=start_hour, minutes=start_minute)
        return date_and_time

    def get_specific_sds_data_set(self, data_set_name):
        return AquaSDSDataSet(self.sd_file_interface.select(data_set_name))

    def get_specific_v_data_set(self, data_set_name):
        return AquaVDataSet(self.v_file_interface.attach(data_set_name))

    def list_data_sets(self):
        for data_set in self.sd_file_interface.datasets().keys():
            print(data_set)

    def get_lat_lon_sets(self):
        latitudes = self.get_specific_sds_data_set("Latitude")
        longitudes = self.get_specific_sds_data_set("Longitude")
        return latitudes, longitudes

    def generate_scans_and_coordinates(self, zones):
        coordinate_list = []
        offnad_scan_list = []
        for area in zones:
            # scan number - 1 = the scan index (from zero)
            # scan index * 2 - the along-track index for MODIS geolocation values (2 values per scan)
            # the +1 ensures the "end" index is at the last of the two geolocation values for the ending scan
            coordinate_list += self.generate_coordinate_data_points((area.get_starting_scan() - 1) * 2,
                                                                    (area.get_ending_scan() - 1) * 2 + 1)
            offnad_scan_list += area.get_scan_list()
        return coordinate_list, offnad_scan_list

    def compare_to_off_nadir(self, nadir_objects, nadir_comp_data, offnad_data="EV_1KM_Emissive"):
        find_boxes_start = time.time()
        print("Finding valid search boxes...")
        search_zones = self.find_zones_with_matches(nadir_objects)
        offn_coords, offn_scans = self.generate_scans_and_coordinates(search_zones)
        find_boxes_end = time.time()
        print("Finihsed! Process took " + str(find_boxes_end-find_boxes_start))
        offnad_data_set = self.get_specific_sds_data_set(offnad_data)
        # 8 = MODIS Band 28. Replace with Index for appropriate band if needed.
        if offnad_data == "EV_1KM_RefSB":
            scale, offset = offnad_data_set.get_scales_and_offsets_for_band(8, "Reflectance")
        else:
            scale, offset = offnad_data_set.get_scales_and_offsets_for_band(8, "Radiance")
        times = self.get_times_list()
        matches = []
        print("Running comparisons...")
        comparison_start = time.time()
        # n - nadir scan index
        # o - off-nadir scan index
        # c - coordinate index (along frame index)
        for n in range(len(nadir_objects)):
            for o in range(len(offn_scans)):
                # offnadir scan time = times[offn_scans[o] - 1]
                if nadir_objects[n].within_time_range(times[offn_scans[o] - 1]):
                    for c in range(len(offn_coords[o])):
                        if nadir_objects[n].within_geospatial_range(offn_coords[o][c]):
                            # c * 5 + 2 - converts the Along Frame Index for a geolocation file into an Along Frame Index for data.
                            tpc = TwoPointComparison(n + 1,
                                                     offn_scans[o],
                                                     nadir_objects[n].get_nadir_pos(),
                                                     (c * 5 + 2),
                                                     offn_coords[o][c],
                                                     nadir_objects[n].get_coordinates())
                            offnad_data_set.compare_values(tpc,
                                                          nadir_comp_data[n],
                                                          s0=scale,
                                                          s1=offset)
                            matches.append(tpc)
        comparison_end = time.time()
        print("Finished! Process took " + str(comparison_end-comparison_start))
        return matches

    def find_zones_with_matches(self, nadir_object_list):
        potential_match_boxes = self.boxes
        times = self.get_times_list()
        valid_zones = []
        for box in potential_match_boxes:
            start_time = times[box.get_starting_scan() - 1]
            end_time = times[box.get_ending_scan() - 1]
            for nadir_point in nadir_object_list:
                if (nadir_point.within_time_range(start_time) or nadir_point.within_time_range(end_time)) and box.encapsulates(
                        nadir_point.get_coordinates()):
                    valid_zones.append(box)
                    break
        return valid_zones

    def find_valid_factor(self):
        # Assumes scans will always range between 202-204 per granule
        scans_factors_dict = {203: 58, 202: 101, 204: 68}
        scans = self.get_number_of_scans()
        return scans_factors_dict[scans]

    def generate_lat_lon_boxes(self):
        latitudes, longitudes = self.get_lat_lon_sets()
        dimensions = latitudes.get_dimensions()
        scale_factor = self.get_scan_to_node_scale_factor(dimensions[0])
        fill_value = latitudes.get_fill_value()
        boxes = []
        div_factor = self.find_valid_factor()
        for i in range(0, dimensions[0], div_factor):
            start_offset = 0
            end_offset = 0
            top_left = (latitudes.get_specific_data_point(i, 0), longitudes.get_specific_data_point(i, 0))
            while top_left == (fill_value, fill_value):
                start_offset += 1
                top_left = (latitudes.get_specific_data_point(i+start_offset, 0), longitudes.get_specific_data_point(i+start_offset, 0))
            top_right = (latitudes.get_specific_data_point(i, dimensions[1]-1), longitudes.get_specific_data_point(i, dimensions[1]-1))
            while top_right == (fill_value, fill_value):
                start_offset += 1
                top_right = (latitudes.get_specific_data_point(i+start_offset, dimensions[1]-1), longitudes.get_specific_data_point(i+start_offset, dimensions[1]-1))
            bottom_right = (latitudes.get_specific_data_point(i + (div_factor-1), dimensions[1]-1), longitudes.get_specific_data_point(i+(div_factor-1), dimensions[1]-1))
            while bottom_right == (fill_value, fill_value):
                end_offset -= 1
                bottom_right = (latitudes.get_specific_data_point(i + (div_factor-1)+end_offset, dimensions[1] - 1), longitudes.get_specific_data_point(i + (div_factor-1)+end_offset, dimensions[1] - 1))
            bottom_left = (latitudes.get_specific_data_point(i + (div_factor-1), 0), longitudes.get_specific_data_point(i + (div_factor-1), 0))
            while bottom_left == (fill_value, fill_value):
                end_offset -= 1
                bottom_left = (latitudes.get_specific_data_point(i + (div_factor-1)+end_offset, 0), longitudes.get_specific_data_point(i + (div_factor-1)+end_offset, 0))
            boxes.append(GeospatialScanBox(top_left, bottom_left, top_right, bottom_right, i // scale_factor + 1, i // scale_factor + 29))
        return boxes

    def generate_coordinate_data_points(self, start_x, end_x):
        latitudes, longitudes = self.get_lat_lon_sets()
        dimensions = latitudes.get_dimensions()
        scale_factor = self.get_scan_to_node_scale_factor(dimensions[0])
        # dimensions[1] - 1 = max y coordinate
        lat_coords = latitudes.get_data_chunk_2d(start_x, end_x, 0, dimensions[1] - 1, scale_factor)
        long_coords = longitudes.get_data_chunk_2d(start_x, end_x, 0, dimensions[1] - 1, scale_factor)
        coordinates = []
        for i in range(len(long_coords)):
            coordinates.append(list(zip(lat_coords[i], long_coords[i])))
        return coordinates

    def get_scan_to_node_scale_factor(self, scaled_dimension):
        number_of_scans = self.get_number_of_scans()
        return scaled_dimension // number_of_scans

    def generate_nadir_data_points(self):
        scan_start_time = self.get_start_time()
        metadata = self.get_specific_v_data_set('Level 1B Swath Metadata')
        # allow inputs for ranges
        nadir_point_list = metadata.generate_nadir_point_search_boxes(scan_start_time, datetime.timedelta(minutes=15), .10)
        return nadir_point_list

    def get_nadir_radiances(self):
        radiances = self.get_specific_sds_data_set('EV_1KM_Emissive')
        # 8 - MODIS Band 28, replace with whichever band.
        ret_vals = radiances.get_nadir_data_by_scan(8)
        return ret_vals

    def get_number_of_scans(self):
        scans = self.attributes['Number of Scans']
        return scans

    def __str__(self):
        return self.name

    def close_file(self):
        self.sd_file_interface.end()
        self.v_file_interface.end()
        self.hdf_file.close()


class HDF5File(object):
    def __init__(self, file_name):
        self.hdf_file = h5py.File(file_name, "r")
        self.main_group = self.hdf_file['All_Data']
        for item in self.main_group.keys():
            if "SDR" in item:
                self.sdr_group = self.main_group[item]
                # 6 = position of band identifier (either I or M) in the SDR folder file name.
                self.file_type = item[6]
            if "GEO" in item:
                self.geo_group = self.main_group[item]
        self.boxes = self.generate_lat_lon_boxes()
        if "RadianceFactors" in self.sdr_group:
            self.c0, self.c1 = self.get_radiance_factors()
        if "ReflectanceFactors" in self.sdr_group:
            self.s0, self.s1 = self.get_reflectance_factors()
        self.name = file_name.split("/")[-1]

    def get_specific_sdr_data_set(self, data_set_name):
        return SuomiDataSet(self.sdr_group[data_set_name])

    def get_specific_geo_data_set(self, data_set_name):
        return SuomiDataSet(self.geo_group[data_set_name])

    def list_data_sets(self):
        for data_set in self.hdf_file.keys():
            print(data_set)

    def get_times_list(self):
        mid_time = self.get_specific_geo_data_set('MidTime')
        return mid_time.convert_all_microseconds()

    def get_lat_lon_sets(self):
        lat = self.get_specific_geo_data_set('Latitude')
        lon = self.get_specific_geo_data_set('Longitude')
        return lat, lon

    def get_number_of_scans(self):
        num_of_scans_set = self.get_specific_geo_data_set('NumberOfScans')
        return num_of_scans_set.sum_single_column_set()

    def get_reflectance_factors(self):
        try:
            rad_factors = self.get_specific_sdr_data_set('ReflectanceFactors')
            c0, c1 = rad_factors.get_scale_factors()
            return c0, c1
        # Default values returned if no data found
        except KeyError:
            return [1,1,1,1], [0,0,0,0]

    def get_radiance_factors(self):
        try:
            rad_factors = self.get_specific_sdr_data_set('RadianceFactors')
            c0, c1 = rad_factors.get_scale_factors()
            return c0, c1
        # Default values returned if no data found
        except KeyError:
            return [1,1,1,1], [0,0,0,0]

    # For consistency, the "nadir along frame index" for the next two functions is the value of the dimension (max value) divded by two.
    def generate_nadir_data_points(self):
        coordinates = self.generate_nadir_coordinates()
        times = self.get_times_list()
        list_of_objs = []
        nadir_along_frame_index = self.get_specific_geo_data_set("Latitude").dimensions[1] // 2
        if len(times) == len(coordinates):
            for scan_no in range(len(coordinates)):
                list_of_objs.append(NadirPoint(coordinates[scan_no][0][0], coordinates[scan_no][1][0], times[scan_no], scan_no,
                                               .10, datetime.timedelta(minutes=15), nadir_along_frame_index))
        return list_of_objs

    def generate_nadir_coordinates(self):
        lat_set, long_set = self.get_lat_lon_sets()
        long_dimensions = long_set.get_dimensions()
        long_coords = long_set.chunk_and_return_scan_data_for(0, long_dimensions[0] - 1, long_dimensions[1] // 2,
                                                              long_dimensions[1] // 2)
        lat_dimensions = lat_set.get_dimensions()
        lat_coords = lat_set.chunk_and_return_scan_data_for(0, lat_dimensions[0] - 1, lat_dimensions[1] // 2,
                                                            lat_dimensions[1] // 2)
        return list(zip(lat_coords, long_coords))

    def compare_to_off_nadir(self, nadir_points, nadir_comp_data, offnad_data="Radiance"):
        print("Finding valid search zones...")
        find_boxes_start = time.time()
        zones = self.find_zones_with_matches(nadir_points)
        offn_coords, offn_scans = self.generate_scans_and_coordinates(zones)
        scan_times = self.get_times_list()
        find_boxes_end = time.time()
        print("Completed in " + str(find_boxes_end - find_boxes_start) + " seconds")
        comparison_set = self.get_specific_sdr_data_set(offnad_data)
        matches = []
        print("Running comparisons ...")
        compare_time_start = time.time()
        for n in range(len(nadir_points)):
            for o in range(len(offn_scans)):
                if nadir_points[n].within_time_range(scan_times[offn_scans[o]-1]):
                    for c in range(len(offn_coords[o])):
                        if nadir_points[n].within_geospatial_range(offn_coords[o][c]):
                            # c*5 is a modification when only every FIFTH element is chosen
                            tpc = TwoPointComparison(offn_scans[o],
                                                     n + 1, c * 5 + 1,
                                                     nadir_points[n].get_nadir_pos(),
                                                     nadir_points[n].get_coordinates(),
                                                     offn_coords[o][c])
                            # while all the scale factors are normally idenitical, the caluclations here ensure that the scale factors for the exact granule are beign used.
                            comparison_set.compare_values(tpc, nadir_comp_data[n],
                                                          self.c0[((offn_scans[o]-1)//48)],
                                                          self.c1[((offn_scans[o]-1)//48)])
                            matches.append(tpc)
        comapre_time_end = time.time()
        print("Completed in " + str(comapre_time_end-compare_time_start) + " seconds")
        return matches

    def generate_scans_and_coordinates(self, geo_zones):
        coordinate_list = []
        scan_list = []
        number_of_detectors = 16
        if self.file_type and self.file_type == "I":
            number_of_detectors = 32
        for area in geo_zones:
            coordinate_list += self.generate_coordinate_data_points(
                (area.get_starting_scan() - 1) * number_of_detectors,
                (area.get_ending_scan() - 1) * number_of_detectors + (number_of_detectors-1))
            scan_list += area.get_scan_list()
        return coordinate_list, scan_list

    def find_zones_with_matches(self, nadir_object_list):
        boxes = self.boxes
        scan_times = self.get_times_list()
        valid_zones = []
        for box in boxes:
            start_scan = box.get_starting_scan()
            end_scan = box.get_ending_scan()
            for item in nadir_object_list:
                if (item.within_time_range(scan_times[start_scan-1]) or item.within_time_range(scan_times[end_scan-1])) and box.encapsulates(item.get_coordinates()):
                    valid_zones.append(box)
                    break
        return valid_zones

    def generate_lat_lon_boxes(self):
        lat_set, long_set = self.get_lat_lon_sets()
        lat_set_max = lat_set.get_dimensions()[1]
        long_set_max = long_set.get_dimensions()[1]
        scan_num = self.get_number_of_scans()
        boxes = []
        for val_index in range(0, scan_num*16, 384):
            start_offset = 0
            end_offset = 0
            top_left = (lat_set.get_specific_data_point(val_index, 0),long_set.get_specific_data_point(val_index, 0))
            while self.is_filler_coordiante(top_left):
                start_offset += 1
                top_left = (lat_set.get_specific_data_point(val_index+start_offset, 0),long_set.get_specific_data_point(val_index+start_offset, 0))
            top_right = (lat_set.get_specific_data_point(val_index, lat_set_max-1), long_set.get_specific_data_point(val_index, long_set_max-1))
            while self.is_filler_coordiante(top_right):
                start_offset += 1
                top_right = (lat_set.get_specific_data_point(val_index+start_offset, lat_set_max-1), long_set.get_specific_data_point(val_index+start_offset, long_set_max-1))
            bottom_right = (lat_set.get_specific_data_point(val_index+383, lat_set_max-1),long_set.get_specific_data_point(val_index+383, long_set_max-1))
            while self.is_filler_coordiante(bottom_right):
                end_offset -= 1
                bottom_right = (lat_set.get_specific_data_point(val_index+383+end_offset, lat_set_max-1),long_set.get_specific_data_point(val_index+383+end_offset, long_set_max-1))
            bottom_left = (lat_set.get_specific_data_point(val_index+383, 0), long_set.get_specific_data_point(val_index+383, 0))
            while self.is_filler_coordiante(bottom_left):
                end_offset -= 1
                bottom_left = (lat_set.get_specific_data_point(val_index+383+end_offset, 0),long_set.get_specific_data_point(val_index+383+end_offset, 0))
            boxes.append(GeospatialScanBox(top_left, bottom_left, top_right, bottom_right, val_index // 16 + 1 + (start_offset//16), val_index // 16 + 24 + (end_offset//16)))
        return boxes

    def is_filler_coordiante(self, coord):
        if coord == (numpy.float32(-999.29999), numpy.float32(-999.29999)):
            return True
        else:
            return False

    def get_nadir_radiances(self):
        scales, offsets = self.get_radiance_factors()
        radiances = self.get_specific_sdr_data_set("Radiance")
        data = radiances.get_nadir_data_by_scan(scales, offsets)
        return data

    def generate_coordinate_data_points(self, start_x, end_x):
        lat_set, long_set = self.get_lat_lon_sets()
        lat_dimensions = lat_set.get_dimensions()
        long_dimensions = long_set.get_dimensions()
        long_coords = long_set.chunk_and_return_scan_data_for(start_x, end_x, 0, long_dimensions[1] - 1)
        lat_coords = lat_set.chunk_and_return_scan_data_for(start_x, end_x, 0, lat_dimensions[1] - 1)
        coordinates = []
        for i in range(len(long_coords)):
            coordinates.append(list(zip(lat_coords[i], long_coords[i])))
        return coordinates

    def __str__(self):
        return self.name

    def close_file(self):
        self.hdf_file.close()