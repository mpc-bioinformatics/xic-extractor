#!/usr/bin/env python

import argparse
import numpy as np
import math
from collections import defaultdict
import datetime
import time
import h5py
import tqdm
import array
import csv
import alphatims
import alphatims.bruker
import bisect
from fisher_py import RawFile
import pandas as pd 

from fisher_py.raw_file_reader import RawFileReaderAdapter, RawFileAccess
from fisher_py.data.business import GenericDataTypes, ChromatogramTraceSettings, TraceType, ChromatogramSignal, SpectrumPacketType, Scan, SegmentedScan, ScanStatistics
from fisher_py.data.filter_enums import MsOrderType
from fisher_py.data import Device, ToleranceUnits
from fisher_py.mass_precision_estimator import PrecisionEstimate


from fisher_py.net_wrapping import ThermoFisher


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-raw", help="RAW Thermo File where to extract XICs")
    parser.add_argument("-query_csv", help="Queries (actual xics) to be retrieved")
    parser.add_argument("-out_hdf5", help="The Output HDF5, containing all XICs")

    return parser.parse_args()

### Modified Python bisect implementations to work with dateframes
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
        if a["rt_values"].iloc[mid] < x: lo = mid+1
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
        if x < a["rt_values"].iloc[mid]: hi = mid
        else: lo = mid+1
    return lo


if __name__ == "__main__":
    args = argparse_setup()

    args.raw = "/home/luxii/Desktop/temp_Svitlana/TIM00011105_S1-F5_1_1341.d"
    args.query_csv = "/home/luxii/Desktop/temp_Svitlana/ISA_SPIKEINS.csv"
    args.out_hdf5 = "/home/luxii/Desktop/temp_Svitlana/delme.hdf5"
    data = alphatims.bruker.TimsTOF(args.raw)

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
        ms1_rt_times = pd.DataFrame(data.convert_from_indices(
                    data[:, :, 0, :, :, "raw"],
                    return_rt_values = True,
                    raw_indices_sorted = True,
                )).groupby(by="rt_values").first().rename_axis('rt_values').reset_index()
        ms2_rt_times = pd.DataFrame({"rt_values": data.rt_values})
        ms2_rt_times = ms2_rt_times[~ms2_rt_times["rt_values"].isin(ms1_rt_times["rt_values"])]
       
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
                rt_start_val = 0
                rt_end_val = data.rt_max_value
            mz_start_val, mz_end_val = float(mz_start_val), float(mz_end_val)
            rt_start_val, rt_end_val = float(rt_start_val)*60, float(rt_end_val)*60  # In Bruker it is in seconds

            # Extract XIC directly
            if ms_val == "ms":
                xic_values = data[rt_start_val:rt_end_val, :, 0, mz_start_val:mz_end_val, :][["rt_values", "corrected_intensity_values"]] \
                    .groupby(by="rt_values").sum().reset_index()

                # Retrieve the scan idcs (set rt window)
                l_scan_idx = bisect_left_rt(ms1_rt_times, rt_start_val, lo=0, hi=len(ms1_rt_times))
                r_scan_idx = bisect_right_rt(ms1_rt_times, rt_end_val, lo=l_scan_idx, hi=len(ms1_rt_times))

                xic_data = ms1_rt_times.iloc[l_scan_idx:r_scan_idx]
                xic_data.loc[:,"intensities"] = 0
                xic_data.loc[xic_data["rt_values"].isin(xic_values["rt_values"]), "intensities"] = xic_values["corrected_intensity_values"].values

            elif ms_val == "ms2":
                xic_values = data[rt_start_val:rt_end_val, :, 1:, mz_start_val:mz_end_val, :][["rt_values", "corrected_intensity_values"]] \
                    .groupby(by="rt_values").sum().reset_index()

                # Retrieve the scan idcs (set rt window)
                l_scan_idx = bisect_left_rt(ms2_rt_times, rt_start_val, lo=0, hi=len(ms1_rt_times))
                r_scan_idx = bisect_right_rt(ms2_rt_times, rt_end_val, lo=l_scan_idx, hi=len(ms1_rt_times))

                xic_data = ms2_rt_times.iloc[l_scan_idx:r_scan_idx]
                xic_data["intensities"] = 0
                xic_data.loc[xic_data["rt_values"].isin(xic_values["rt_values"]), "intensities"] = xic_values["corrected_intensity_values"].values

            else:
                raise Exception("Cannot extract ms level {}".format(ms_val))

            # Save in h5
            out_h5["retention_times"][h5_idx] = array.array("d", xic_data["rt_values"])
            out_h5["intensities"][h5_idx] = array.array("d", xic_data["intensities"])
