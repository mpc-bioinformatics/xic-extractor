#!/usr/bin/env python

import argparse
import array
import bisect
import csv

import h5py
import numpy as np
import tqdm
from fisher_py.data import Device
from fisher_py.raw_file_reader import RawFileReaderAdapter


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-raw", help="RAW Thermo File where to extract XICs")
    parser.add_argument("-query_csv", help="Queries (actual xics) to be retrieved")
    parser.add_argument("-out_hdf5", help="The Output HDF5, containing all XICs")

    return parser.parse_args()

### Modified Python bisect implementations to work with 2 dimensional arrays
def bisect_left_rt(a, x, lo=0, hi=None):
    """Return the index where to insert item x in list a, assuming a is sorted.

    The return value i is such that all e in a[:i] have e < x, and all e in
    a[i:] have e >= x.  So if x already appears in the list, a.insert(x) will
    insert just before the leftmost x already there.

    Optional args lo (default 0) and hi (default len(a)) bound the
    slice of a to be searched.
    """
    if lo < 0:
        raise ValueError('lo must be non-negative')
    if hi is None:
        hi = a.shape[0]
    while lo < hi:
        mid = (lo+hi)//2
        if a(mid) < x: lo = mid+1
        else: hi = mid
    return lo


def bisect_right_rt(a, x, lo=0, hi=None):
    """Return the index where to insert item x in list a, assuming a is sorted.

    The return value i is such that all e in a[:i] have e <= x, and all e in
    a[i:] have e > x.  So if x already appears in the list, a.insert(x) will
    insert just after the rightmost x already there.

    Optional args lo (default 0) and hi (default len(a)) bound the
    slice of a to be searched.
    """
    if lo < 0:
        raise ValueError('lo must be non-negative')
    if hi is None:
        hi = a.shape[0]
    while lo < hi:
        mid = (lo+hi)//2
        if x < a(mid): hi = mid
        else: lo = mid+1
    return lo


if __name__ == "__main__":
    args = argparse_setup()

    raw_file = RawFileReaderAdapter.file_factory(args.raw)
    raw_file.select_instrument(Device.MS, 1)  # Selecting the MS

    first_scan_number = raw_file.run_header_ex.first_spectrum
    last_scan_number = raw_file.run_header_ex.last_spectrum

    with h5py.File(args.out_hdf5, "w") as out_h5, open(args.query_csv, "r") as q_in:
        # Read CSV - Queries
        q_csv = csv.reader(q_in)
        headers = next(q_csv)

        # Get the corresponding headers
        ident_idx = headers.index("identifier")
        mz_start_idx = headers.index("mz_start")
        mz_end_idx = headers.index("mz_end")
        rt_start_idx = headers.index("rt_start")
        rt_end_idx = headers.index("rt_end")
        mz_idx = headers.index("mz")
        ppm_idx = headers.index("ppm")
        ms_idx = headers.index("ms_level")

        # Load all scans, ms_level and queries into memory
        queries = [l for l in q_csv]  # Process Queries, as they have been added
        ms_level = [x.ms_order.name.lower() for x in raw_file.get_scan_events(first_scan_number, last_scan_number)]
        all_scans_intens = [np.array(raw_file._get_wrapped_object_().GetSegmentedScanFromScanNumber(x, None).Intensities) for x in  range(first_scan_number, last_scan_number)]
        all_scans_pos = [np.array(raw_file._get_wrapped_object_().GetSegmentedScanFromScanNumber(x, None).Positions) for x in  range(first_scan_number, last_scan_number)]

        # Set output with variable length array
        dt = h5py.vlen_dtype(np.dtype('float64'))
        out_h5.create_dataset("retention_times", (len(queries),), dtype=dt, compression="gzip")
        out_h5.create_dataset("intensities", (len(queries),), dtype=dt, compression="gzip")
        out_h5.create_dataset("labels", data=[x[ident_idx] for x in queries], compression="gzip")

        # Retrieve for each single query:
        for h5_idx, l in enumerate(tqdm.tqdm(queries, unit="queries")):
            # Initialize values
            ident_val = l[ident_idx]
            mz_start_val = l[mz_start_idx]
            mz_end_val = l[mz_end_idx]
            rt_start_val = l[rt_start_idx]
            rt_end_val = l[rt_end_idx]
            mz_val = l[mz_idx]
            ppm_val = l[ppm_idx]
            ms_val = l[ms_idx].replace("/", "")
            if not mz_start_val:
                # We retrieve mz start and end via ppm error tolerance
                tol = float(mz_val) / 1000000 * float(ppm_val)
                mz_start_val = float(mz_val) - tol
                mz_end_val = float(mz_val) + tol
            if not rt_start_val:
                # If no retention time is given, return for all rts
                rt_start_val = raw_file.run_header_ex.start_time
                rt_end_val = raw_file.run_header_ex.end_time
            mz_start_val, mz_end_val = float(mz_start_val), float(mz_end_val)
            rt_start_val, rt_end_val = float(rt_start_val), float(rt_end_val)

            # Retrieve the scan idcs (filter rt)
            l_scan_idx = bisect_left_rt(raw_file.retention_time_from_scan_number, rt_start_val, lo=first_scan_number, hi=last_scan_number)
            r_scan_idx = bisect_right_rt(raw_file.retention_time_from_scan_number, rt_end_val, lo=l_scan_idx, hi=last_scan_number)

            # Retrieve the XIC (filter mz)
            xic = np.zeros((r_scan_idx - l_scan_idx, 2))
            bool_arr = np.zeros((r_scan_idx - l_scan_idx,), dtype=bool)
            for idx, i in enumerate(range(l_scan_idx, r_scan_idx)):
                if ms_level[i-1] == ms_val:
                    scan_pos = all_scans_pos[i-1]
                    scan_intens = all_scans_intens[i-1]

                    l_pos = bisect.bisect_left(scan_pos, mz_start_val, lo=0, hi=scan_pos.shape[0])
                    r_pos = bisect.bisect_right(scan_pos, mz_end_val, lo=l_pos, hi=scan_pos.shape[0])
                    xic[idx, 0] = raw_file.retention_time_from_scan_number(i)
                    xic[idx, 1] = scan_intens[l_pos:r_pos].max(initial=0)  # Get the BasePeak if multiple are present!
                    bool_arr[idx] = True

            # Filter values out for skipped entries
            xic = xic[bool_arr]

            # Save in h5
            out_h5["retention_times"][h5_idx] = array.array("d", xic[:,0])
            out_h5["intensities"][h5_idx] = array.array("d", xic[:,1])
