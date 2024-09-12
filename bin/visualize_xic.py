#!/bin/env python

import argparse
import csv
import sys
import json
import os
import pyopenms
import re
import h5py
import pandas as pd
import plotly
import plotly.graph_objects as go
import plotly.express as ex
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-group_regex", help="If wanted, use a regex expression to group the file. Defaults, to just use the whole stirng (each file individually). The first captured group in Regex will be used for grouping", default="(.*)")
    parser.add_argument("-input_hdf5s", help="The json files (in TRFP XIC-Extraction Format), seperated by comma: ASSUMPTION, labels are all sorted in the same order across hdf5s!!")
    parser.add_argument("-input_trafos", help="The corresponding_trafo.xml files (from a Openms RT-Alignment), seperated by comma")
    parser.add_argument("-outdir", help="The output folder in which all the figures should be saved.")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    # Filter and get data lists
    input_hdf5s = args.input_hdf5s.split(",")
    # Extra filter step for the optional trafo xmls
    input_trafos = [x for x in args.input_trafos.split(",") if not x.startswith("NULL___")]

    data_dict = dict()
    for hdf5_file in input_hdf5s:
        filename = ".".join(os.path.basename(hdf5_file).split(".")[:-1])

        o_h5 = h5py.File(hdf5_file, 'r')

        # Get Grouping info (if any)
        matches = re.finditer(args.group_regex, filename, re.MULTILINE)
        for matchNum, match in enumerate(matches, start=1):
            group = (match.group(0))
            print(group)
            break

        # Content:
        # dict[FILENAME] --> [GROUP, TimeMultiplier, h5file, OPTIONAL[trafo.xml]]
        data_dict[filename] = [
            group,
            60,  # Time unit is in minutes, we want to visualize seconds
            o_h5,
            None
            ]

    # If available: read the Trafo
    for trafo_file in input_trafos:
        filename = ".".join(os.path.basename(trafo_file).split(".")[:-1])
        if filename in data_dict:
            # Read trafo file via Openms
            trafo_xml = pyopenms.TransformationXMLFile()
            trafo_desc = pyopenms.TransformationDescription()
            trafo_xml.load(trafo_file, trafo_desc, True)

            # Set it in the dict
            data_dict[filename][3] = trafo_desc

    # Get all extracted xic labels
    all_labels = set()
    for key, val in data_dict.items():
        all_labels = all_labels.union(
            set(val[2]["labels"])
        )

    for idx, label in enumerate(data_dict[next(iter(data_dict))][2]["labels"]):
        identifier = label.decode("utf-8")

        ### Visualize plots
        # First create a Pandas Dataframe: 
        df_list = []  # [(FileName, Group, RT, Intensity)]
        for key, value in data_dict.items():
            filename = [key]*len(value[2]["retention_times"][idx])
            group = [value[0]]*len(value[2]["retention_times"][idx])
            if value[-1] is not None:
                # Apply RT Transformation
                rts = [value[-1].apply(x*value[1]) for x in value[2]["retention_times"][idx]]
            else:
                # Just copy the original values    
                rts = [x*value[1] for x in value[2]["retention_times"][idx]]
            intens = value[2]["intensities"][idx]
            df_list += [(w,x,y,z) for w,x,y,z in zip(filename, group, rts, intens)]

        df = pd.DataFrame(df_list, columns=["File", "Group", "RT", "Intens"])


        # Plot Information in one plot
        color_map = {y: ex.colors.qualitative.Plotly[x] for x, y in enumerate(set(df["Group"]))}
        fig = go.Figure()
        for key in data_dict.keys():
            fig.add_trace(go.Scatter(
                x=df[df["File"] == key]["RT"],
                y=df[df["File"] == key]["Intens"],
                legendgroup=df[df["File"] == key]["Group"].iloc[0],
                legendgrouptitle_text="Group: " + df[df["File"] == key]["Group"].iloc[0],
                name=key,
                mode="lines",
                line=dict(color=color_map[df[df["File"] == key]["Group"].iloc[0]])
            ))
        fig.update_layout(title=identifier)
        fig.update_layout(legend=dict(groupclick="toggleitem"))

        plotly.offline.plot(fig, filename=args.outdir + os.sep + identifier + ".html", auto_open=False)
        fig.write_image(args.outdir + os.sep + identifier + ".png")
