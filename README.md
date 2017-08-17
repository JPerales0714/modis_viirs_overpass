# VIIRS / MODIS Overpass Analysis Tool

Primary purpose of project is to provide a framework for identifying, collecting, and analyzing data at (near) simultaneous overpasses between Aqua (MODIS) and Suomi-NPP VIIRS.

Most tests were run using MODIS Band 28 and VIIRS Band M14 Radiance values. These are the defaults in the code presently.

Additional functionality provided for creating heatmaps that highlight differences between values picked up by VIIRS I-Bands and VIIRS M-Bands (I1 and M5 specifically).

## Getting Started

### Prerequisites

To run this code, you, need to install:

* Python 3+ (Required)
* Anaconda (Highly recommended for installing the other extensions)
* PyHDF (Required)
* H5PY (Required) 
* Psycopg2 (If you use a PostgreSQL database)
* PyCharm (Recommended)

### Data Repositories Used

I highly recommend using [LAADS DAAC] (https://ladsweb.modaps.eosdis.nasa.gov/) for MODIS data and [CLASS] (https://www.class.ncdc.noaa.gov/saa/products/welcome) for VIIRS data. Code was written with the file formats available from these sites in mind.

## Running the program

Assuming your computer is configured correctly with all the correct libraries, you only need to do take note of the following:

* The "Base File Path" should be an absolute path to your data directory.
* All files to be used must be in a seperate directory (folder), the code doesn't support singular files as inputs.
* Tables within a database are not generated as part of this code. Any specified tables to be used must be pre-existig, or table initialization code must be added.
* Similarly, data tables are not truncated unless commands are added.
* HDF4 and HDF5 files are required, but either type can be used as Nadir or Off-Nadir data.

Otherwise, simply follow the prompting instructions on-screen.

## Notes regarding design

* Changing MODIS Bands to analyze or the temporal/spatial search criteria for a match must be modified within the code.
* Several functions are abstracted to support several types of inputs, but the program wasn't necessarily designed to support the same, so be cautious of results from anything save Longwave to Longwave band comparisons. Existing functions can always be utilized in different ways.
* While Frame Positions (referred to as Along Track Indices) are always in their "from zero" (indexable) format, Scans are frequently in their "numerical" (counted) format. As a result, whenever indexing using scans, 1 must be subtracted from the scan value to become the correct corresponding index. The benefit of this is that scan values can be printed and easily understood. 
* Currently info_to_database averages all variances (of radiance values) per scan angle and submits those points to the database. This function can be modified to utilize any and all of the data from the Two-Point-Comparison objects.
* This code was meant to be introductory - so certain error-handling operations, opprotunities for shorter code, and frivilous method defining was ignored.

## Contact Info

For any questions, please contact jmp3zb@virginia.edu
