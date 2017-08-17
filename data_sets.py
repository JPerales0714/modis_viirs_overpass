from pyhdf.SD import SDS
from pyhdf.VS import *
from pyhdf.HDF import *
import datetime
import numpy
from data_structures import NadirPoint


class AquaVDataSet(object):
    def __init__(self, data_set):
        if isinstance(data_set, VD):
            self.data = data_set
        else:
            raise Exception("Invalid data set type for AquaVDataSet")

    def generate_nadir_point_search_boxes(self, scan_start_time, search_box_duration, search_box_space):
        point_list = []
        scan_time = scan_start_time
        while 1:
            try:
                record = self.data.read()
                latitude = record[0][7]
                longitude = record[0][8]
                scan_number = record[0][0]
                nadir_along_swath_frame = record[0][6]
                nadir_point = NadirPoint(latitude, longitude, scan_time, scan_number, search_box_space,
                                         search_box_duration, nadir_along_swath_frame)
                point_list.append(nadir_point)
                scan_time += datetime.timedelta(seconds=1.4771)
            except HDF4Error:
                break
        return point_list


class AquaSDSDataSet(object):

    def __init__(self, data_set):
        if isinstance(data_set, SDS):
            self.data = data_set
            self.attributes = data_set.attributes()
            self.info = data_set.info()
            self.rank = self.info[1]
            self.name = self.info[0]
            if "KM" in self.info[0]:
                # dictionary describing types of MODIS files -> number of detectors
                detectors = {"1KM": 10, "QKM": 40, "HKM": 20}
                name_list = self.info[0].split("_")
                self.num_of_detectors = detectors[name_list[1]]
        else:
            raise Exception("Invalid data set type for AquaSDSDataSet")

    def get_attributes(self):
        return self.attributes

    def get_dimensions(self):
        return self.info[2]

    def get_fill_value(self):
        return self.attributes["_FillValue"]

    def get_specific_data_point(self, x, y):
        return float(self.data.get([x, y], [1, 1]))

    def compare_values(self, match, viirs_value, s0=1, s1=1):
        scan_value = match.get_modis_scan()
        swath_pos = match.get_modis_swath_pos()
        # Scan value is NOT scaled from zero, so one must be subtracted.
        modis_base = self.get_data_chunk_3d(8,
                                            (scan_value - 1) * self.num_of_detectors,
                                            scan_value * self.num_of_detectors - 1,
                                            swath_pos,
                                            swath_pos,
                                            self.num_of_detectors)
        modis_adjusted = s0 * (modis_base[0][0].item() - s1)
        match.set_comparison_values_modis_offnad(viirs_value, modis_adjusted, swath_pos)

    def get_data_chunk_3d(self, band, start_x, end_x, start_y, end_y, number_of_values_per_scan):
        if self.rank == 3:
            final_list = []
            big_list = []
            data_subset = self.data[band]
            filter_subset = data_subset[start_x:end_x + 1, start_y:end_y + 1]
            for row in filter_subset:
                big_list.append(row)
                if len(big_list) % number_of_values_per_scan == 0:
                    final_list.append((list(map(self.data_mean, zip(*big_list)))))
                    big_list = []
            return final_list
        else:
            raise Exception("Attempted to 3D-Chunk a Non-3D data set.")

    def get_data_chunk_2d(self, start_x, end_x, start_y, end_y, number_of_values_per_scan):
        if self.rank == 2:
            final_list = []
            big_list = []
            data_subset = self.data.get([start_x, start_y], [(end_x - start_x) + 1, (end_y - start_y) + 1])
            for row in data_subset:
                big_list.append(row)
                if len(big_list) % number_of_values_per_scan == 0:
                    # zip(*big_list) takes all the lists (rows) in big_list and creates one new list, where each element ...
                    # is another list containing the values from each row at corresponding indices.
                    # [1,2,3] and [2,3,4] zipped becomes [[1,2],[2,3],[3,4]]
                    # map applies the 'data_mean' function to every member of zip(*big_list)
                    final_list.append(list(map(self.data_mean, zip(*big_list))))
                    big_list = []
            return final_list
        else:
            raise Exception("Attempted to 2D-Chunk a Non-2D data set.")

    def data_mean(self, a):
        return sum(a) / len(a)

    def get_scales_and_offsets_for_band(self, band, dtype):
        if dtype == "Radiance":
            scales = self.attributes['radiance_scales']
            offsets = self.attributes['radiance_offsets']
            return float(scales[band]), float(offsets[band])
        elif dtype == "Reflectance":
            scales = self.attributes['reflectance_scales']
            offsets = self.attributes['reflectance_offsets']
            return float(scales[band]), float(offsets[band])

    def get_nadir_data_by_scan(self, band):
        nadir_frame = self.info[2][2] // 2
        along_track_len = self.info[2][1]
        if "Emissive" in self.name:
            scale, offset = self.get_scales_and_offsets_for_band(band, "Radiance")
        elif "RefSB" in self.name:
            scale, offset = self.get_scales_and_offsets_for_band(band, "Reflectance")
        else:
            # defaults
            scale, offset = 1, 0
        data = {}
        for track_number in range(0, along_track_len, self.num_of_detectors):
            # Band/Track/Frame
            row = self.data.get([band, track_number, nadir_frame], [1, self.num_of_detectors, 1])
            total = 0
            for i in row.flat:
                total += i
            avg = total // self.num_of_detectors
            adjusted_val = scale * (avg - offset)
            data[track_number // self.num_of_detectors] = adjusted_val
        return data


class SuomiDataSet(object):

    def __init__(self, data_set):
        if isinstance(data_set, SDS):
            raise Exception("Invalid data set type for SuomiDataSet")
        else:
            self.ref_data = data_set
            self.data = numpy.array(self.ref_data)
            self.attributes = self.ref_data.attrs
            self.dimensions = self.ref_data.shape
            # the -1 works on even 1D data sets
            if self.dimensions[-1] % 3200 == 0:
                if self.dimensions[1] // 3200 == 1:
                    self.num_of_detectors = 16
                    self.band_type = "M"
                if self.dimensions[1] // 3200 == 2:
                    self.num_of_detectors = 32
                    self.band_type = "I"

    def get_dimensions(self):
        return tuple(self.dimensions)

    def get_attributes(self):
        return self.attributes

    def get_specific_data_point(self, x, y):
        return self.data[x, y]

    def get_scale_factors(self):
        c0 = []
        c1 = []
        for i in range(0, len(self.data), 2):
            c0.append(self.data[i])
            c1.append(self.data[i + 1])
        return c0, c1

    def chunk_and_return_scan_data_for(self, start_x, end_x, start_y, end_y):
        final_list = []
        big_list = []
        data_subset = self.data[start_x:(end_x + 1), start_y:(end_y + 1)]
        for row in data_subset:
            big_list.append(self.get_elements_at_interval(row, 5))
            if len(big_list) % self.num_of_detectors == 0:
                final_list.append(list(map(self.data_mean, zip(*big_list))))
                big_list = []
        return final_list

    # This function is for I-Band data sets (Reflectance, Radiances) ONLY.
    def get_aggregate_value(self, ref_x, ref_y):
        if self.band_type == "I":
            data_subset = self.data[ref_x:ref_x + 2, ref_y:ref_y + 2]
            sum = 0.0
            for row in data_subset:
                for i in row:
                    sum += i
            return sum / 4.0
        else:
            raise Exception("4x4 Aggregation only intended for I-Band data sets.")

    def get_elements_at_interval(self, unmodified_list, interval):
        modified_list = []
        for i in range(0, len(unmodified_list), interval):
            modified_list.append(unmodified_list[i])
        return modified_list

    def sum_single_column_set(self):
        num = 0
        for column in self.data:
            num += column
        return num

    def data_mean(self, a):
        return sum(a) / len(a)

    def convert_all_microseconds(self):
        # Inteded for: MidTime and StartTime sets.
        times = []
        for time_value in self.data:
            time_value = time_value.item()
            value_in_microseconds = datetime.timedelta(microseconds=time_value)
            corrected_time = datetime.datetime(1958, 1, 1) + value_in_microseconds
            times.append(corrected_time)
        return times

    def compare_values(self, match, modis_value, s0=1, s1=1):
        scan_value = match.get_viirs_scan()
        swath_pos = match.get_viirs_swath_pos()
        # scan value is NOT scaled form zero, so 1 must be subtracted.
        viirs_base = self.chunk_and_return_scan_data_for((scan_value - 1) * self.num_of_detectors,
                                                         scan_value * self.num_of_detectors - 1,
                                                         swath_pos,
                                                         swath_pos)
        viirs_adjusted = s0 * viirs_base[0][0].item() + s1
        match.set_comparison_values_viirs_offnad(viirs_adjusted, modis_value, swath_pos)

    def get_nadir_data_by_scan(self, scales, offsets):
        nadir_frame = self.dimensions[1] // 2
        scans_times16 = self.dimensions[0]
        radiances = {}
        for track_number in range(0, scans_times16, 16):
            scale = scales[track_number // (48 * self.num_of_detectors)]
            offset = offsets[track_number // (48 * self.num_of_detectors)]
            row = self.data[track_number:track_number + self.num_of_detectors, nadir_frame]
            total = 0
            for i in row.flat:
                total += i
            avg = total // self.num_of_detectors
            adjusted_val = scale * avg + offset
            radiances[track_number // self.num_of_detectors] = adjusted_val
        return radiances
