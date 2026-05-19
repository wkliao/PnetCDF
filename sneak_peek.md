------------------------------------------------------------------------------
This is essentially a placeholder for the next release note ...
------------------------------------------------------------------------------

* New feature
  + Intra-node aggregation for read requests is added. This is the complement
    to the write counterpart first implemented in version 1.14.0. Now the
    intra-node aggregation feature supports both write and read operations.
    This feature can be enabled by setting hint `nc_num_aggrs_per_node` to the
    desired number of aggregators per compute node.

* New optimization
  + When creating a new file on the Lustre file system, PnetCDF will try to set
    the Lustre file striping count to the number of compute nodes (NUMA nodes)
    allocated to the MPI job. See the section of new PnetCDF hints below for
    detailed information.

* New I/O driver
  + A new internal I/O driver, named "GIO", is added. It is an alternative I/O
    driver to the MPI-IO driver. GIO includes several implementations for
    performance improvement, such as:
    * automatically sets `cb_nodes` if this hint is not set by the user.
    * supports `cb_nodes` to be set to a multiple of the number of OSTs on
      Lustre.
    * supports hint `cb_buffer_size` to be a multiple of file striping unit
      size.
    * supports Lustre overstriping through hint `overstriping_ratio`.

* New Limitations
  + none

* Configure options

* Configure updates:
  + none

* New constants
  + none

* New APIs
  + none

* APIs deprecated
  + The "vard" APIs introduced in version 1.6.0 are now deprecated. These are
    the API family that take an argument of MPI derived data type describing
    the file access layout, which is used as the fileview by the underlying
    MPI-IO library. The deprecation is because direct file access can be a
    security risk and error prone.

* API syntax changes
  + none

* API semantics updates
  + API `ncmpi_inq_header_size()` now can be called in the define mode. This
    API returns the file header size with metadata defined by the time of the
    call. The inquired file header size can be helpful to pick proper values
    for arguments `h_minfree`, `v_align`, `v_minfree`, `r_align` when calling
    API `ncmpi__enddef()` to make a sufficiently large free space for file
    header extent and variable data sections to grow without moving data
    already stored in the file, i.e. when adding new variables, dimensions, or
    attributes.
    See [PR #201](https://github.com/Parallel-NetCDF/PnetCDF/pull/201).

* New error code precedence
  + none

* Updated error strings
  + none

* New error code
  + `NC_EFSTYPE` indicates an error when an invalid file system type is
    detected by the GIO library.
  + `NC_EDRIVER` indicates an invalid PnetCDF I/O driver is set in I/O hint.
    The current valid drivers are "gio" and "mpiio".
  + `NC_EFILEVIEW` indicates a PnetCDF internal error when creating an MPI
    fileview whose offsets violate the MPI standard requirement of being in a
    monotonically non-decreasing order.

* New PnetCDF hints
  + `nc_data_move_chunk_size` -- When adding new data objects into an existing
    file, the data sections may need to be moved to a higher file offset. The
    movement is performed in chunks. This hint allows users to customized the
    chunk size. The default is 1048576 bytes, i.e. 1 MiB.
    See [PR #203](https://github.com/Parallel-NetCDF/PnetCDF/pull/203).
  + `nc_file_striping` -- When creating a new file on the Lustre file system,
    this hint advises PnetCDF to set the new file's striping configuration. The
    hint value is either "auto" or "inherit". The former sets the new file's
    striping unit to 1 MiB and striping count to the number of compute nodes
    found in the MPI communicator passed to `ncmpi_create()`. The latter sets
    the striping of new file to inherit its parent folder's striping settings,
    if the folder's striping is set. However, if users also set the MPI-IO hint
    `striping_factor` or `striping_unit`, then these MPI-IO hints will take a
    higher precedence. This hint's default value is "auto".
    See [PR #222](https://github.com/Parallel-NetCDF/PnetCDF/pull/222).
  + `nc_driver` -- To select the internal I/O driver. String value of "gio" is
    to use the GIO driver, an external library used as a sub-module of PnetCDF,
    and "mpiio" the MPI-IO driver. The default is "gio".

* New run-time environment variables
  + none

* Build recipes
  + none

* Updated utility programs
  + none

* Other updates:
  + none

* Bug fixes
  + Fix data movement when new record variables are added to an existing file
    that does not change the starting offset of record variable section.
    See [PR #199](https://github.com/Parallel-NetCDF/PnetCDF/pull/199).

* New example programs
  + none

* New programs for I/O benchmarks
  + none

* New test program
  + test/testcases/tst_grow_data.c -- adding new variables by re-entering the
    define mode multiple time, but does not cause file header extent to grow.
    It also tests a case when adding a new record variable that does not change
    the starting offset of the record variable section in the file.

* Issues with NetCDF library
  + none

* Conformity with NetCDF library
  + none

* Discrepancy from NetCDF library
  + none

* Issues related to MPI library vendors:
  + none

* Issues related to Darshan library:
  + none

* Clarifications about of PnetCDF hints
  + none

