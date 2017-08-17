import os
from file_handler import HDF4File, HDF5File
import matplotlib.pyplot as plt
import psycopg2
import numpy
import pickle


base_location_file = 'basepath.pk'
database_info_file = 'dbinfo.pk'
base_location = ''
base_db_info = []

def run_program():
    print("Current Data Analysis Functionalities: ")
    print(" - Nadir vs. Off-Nadir Long-Wave Data Comparison [NVON]")
    print(" - I vs. M-Band Image Overlays for VIIRS [IVM]")
    response = input("Which operations would you like to perform? (NVON/IVM): ")
    if response.lower() == "nvon":
        input_directory_info()
        print("For the nadir data ... ")
        nadir_files = gather_input_files()
        print("For the off-nadir data ...")
        off_nadir_files = gather_input_files()
        input_db_info()
        response2 = input("Would you like to run the reverse as well? [y/n]: ")
        nvon_options_and_run(nadir_files, off_nadir_files)
        # If response2 is a faulty input, the reverse is simply not executed.
        if response2.lower() == "y":
            nvon_options_and_run(off_nadir_files, nadir_files)
    elif response.lower() == "ivm":
        input_directory_info()
        print("Note: for I and M-Band inputs, only the first files in the specified folders will be processed.")
        print("For the I1 data ... ")
        i_files = gather_input_files()
        print("For the M5 data ... ")
        m_files = gather_input_files()
        ivm(m_files[0], i_files[0])
    else:
        print("Invalid entry, please use one of the acronyms listed below:")
        run_program()

def input_directory_info():
    global base_location_file
    global base_location
    resp = input("Would you to change the base directory for inputs? [y/n]: ")
    if resp.lower() == "y":
        location = input("Please input the new path: ")
        while not os.path.isdir(location):
            location = input("Invalid path - please re-enter directory: ")
        base_location = location
        with open(base_location_file, 'wb') as file:
            pickle.dump(base_location, file)
    elif resp.lower() == "n":
        with open(base_location_file, "rb") as file:
            base_location = pickle.load(file)
    else:
        print("Invalid inputs - please input 'y' or 'n'")
        input_directory_info()

def gather_input_files():
    global base_location
    folder = input("Please specify the directory for the data: ")
    path = base_location + "/" + folder
    while not os.path.isdir(path):
        folder = input("Invalid path - please re-enter directory: ")
        path = base_location + "/" + folder
    file_list = []
    for filename in os.listdir(path):
        if filename.endswith(".h5") or filename.endswith("hdf"):
            print("Opening ... " + filename)
            file_list.append(open_file(path + "/" + filename))
    return file_list

def open_file(given_file):
    f_name = given_file
    if f_name[-3:] == "hdf":
        opened_file = HDF4File(f_name)
        return opened_file
    elif f_name[-2:] == "h5":
        opened_file = HDF5File(f_name)
        return opened_file

def input_db_info():
    global base_db_info
    global database_info_file
    resp = input("Would you like to change the database information? [y/n]: ")
    if resp.lower() == "y":
        conn = False
        while conn == False:
            print("Please enter the following: ")
            db_name = input("Database name: ")
            user_name = input("User name: ")
            password = input("Password: ")
            try:
                psycopg2.connect("dbname="+db_name+" user="+user_name+" password="+password)
            except psycopg2.Error as e:
                print("Unable to connect!")
                print("Error: " + str(e.pgerror))
            # This will only be reached once the "try" statement executes successfuly.
            else:
                conn = True
                base_db_info = [db_name, user_name, password]
                with open(database_info_file, 'wb') as file:
                    pickle.dump(base_db_info, file)
    elif resp.lower() == "n":
        with open(database_info_file, "rb") as file:
            base_db_info = pickle.load(file)
    else:
        print("Invalid inputs - please input 'y' or 'n'")
        input_db_info()

def info_to_database(list_of_matches, table, cursor):
    final_dict = {}
    for pair in list_of_matches:
        if .9 <= pair.get_ratio() <= 1.1:
            if pair.get_angle() not in final_dict.keys():
                final_dict[pair.get_angle()] = [pair.get_ratio()]
            else:
                final_dict[pair.get_angle()].append(pair.get_ratio())
    for key in final_dict.keys():
        cursor.execute("INSERT INTO " + table + " (angle, avgv, std) VALUES (%s, %s, %s)",
                    (key, numpy.mean(final_dict[key]), numpy.std(final_dict[key])))

def nvon_options_and_run(nadir_files, on_files):
    global base_db_info
    # base_db_info = [database, username, password]
    db = input("Would you like to submit this data to your database? [y/n]: ")
    # Submitting to database initialized to false.
    database_submit = False
    # Initialization / default for the table name
    table = ""
    if db.lower() == "y":
       database_submit = True
    # Input-checking not done for table names due to limitations in interacting with the database.
       table = input("What table would you like to submit the data to?: ")
    conn = psycopg2.connect("dbname=" + base_db_info[0] + " user=" + base_db_info[1] + " password=" + base_db_info[2])
    cur = conn.cursor()
    final_dict = {}
    for n_num in range(len(nadir_files)):
        n_points = nadir_files[n_num].generate_nadir_data_points()
        n_radiances = nadir_files[n_num].get_nadir_radiances()
        for on_num in range(len(on_files)):
            print("Comparing off-nadir " + str(on_files[on_num]) + " to nadir values of " + str(nadir_files[n_num]))
            final_dict[(on_num, n_num)] = on_files[on_num].compare_to_off_nadir(n_points, n_radiances)
            print("Found " + str(len(final_dict[(on_num, n_num)])) + " matches.")
            if database_submit:
                info_to_database(final_dict[(on_num, n_num)], table, cur)
    conn.commit()
    cur.close()
    conn.close()

def ivm(m_file, i_file):
    trace = []
    m = []
    i = []
    # To avoid constantly re-loading data into memory, data set information is extracted for this function.
    i_data = i_file.get_specific_sdr_data_set("Reflectance")
    m_data = m_file.get_specific_sdr_data_set("Reflectance")
    f1, f2 = i_file.get_reflectance_factors()
    s1, s2 = m_file.get_reflectance_factors()
    # 3072 = 48 (number of scans in a granule) x 4 (number of granules in a CLASS file) x 16 (number of detectors per scan for M Files).
    for x in range(3072):
        row = []
        mrow = []
        irow = []
        # 3200 (standard number of frames for an M file)
        for y in range(3200):
            # Scale Factors f1/f2 and s1/s2 are indexed per granule. x position // number of indexes per granule = granule index
            i_value = (i_data.get_aggregate_value(x * 2, y * 2) * f1[x // 1536] + f2[x // 1536]).item()
            m_value = (m_data.get_specific_data_point(x, y) * s1[x // 768] + s2[x // 768]).item()
            # i1 - .93*m5
            row.append((i_value) - (.93*m_value))
            mrow.append(m_value)
            irow.append(i_value)
        trace.append(row)
        m.append(mrow)
        i.append(irow)
    create_heatmap(trace, "I1-M5")
    create_heatmap(i, "I1")
    create_heatmap(m, "M5")
    plt.show()

def create_heatmap(matrix, title):
    plt.figure()
    plt.title(title)
    plt.imshow(matrix)
    plt.colorbar()

if __name__ == '__main__':
    run_program()