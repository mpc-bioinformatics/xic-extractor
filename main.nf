#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters required for standalone execution
params.raw_spectra = "$PWD/raws"  // RAW-SPectra (either .d (TODO) or .RAW)
params.extraction_csv = "$PWD/xics_to_extract.csv" // A CSV list (as described in the README.md), which tells this worklfow which XICS (or other slicings) to extract
params.group_regex = "(.*)"  // Regex for grouping. This is used for coloring in the visualizations. Leave empty for no coloring
params.export_visualization = false

// Optional Parameters
params.trafoxmls = "" // RT aligned XMLs (from OpenMS) which need to correspond to the raw_spectra files (for RT alignment). Leave empty for alignment
params.outdir = "$PWD/results"  // Output-Directory of the XICs and visualizations


workflow {
    rawfiles = Channel.fromPath(params.raw_spectra + "/*.{raw,d}", type: "any")
    xics_to_extract = Channel.fromPath(params.extraction_csv)
    if(params.trafoxmls == "") {
        trafoXMLs = Channel.empty()
    } else {
        trafoXMLs = Channel.fromPath(params.trafoxmls + "/*.trafoXML")
    }
    
    extract_xics(rawfiles, xics_to_extract, trafoXMLs)
}

workflow extract_xics {
    take:
        raw_files  // RAW spectra files (.raw or .d)
        xics_to_extract  // CSV List of the XICs to be extracted.
        trafoXMLs  // Optional Parameter. Provide "Channel.empty()" if no trafoXMLs are present
    main:
        // Retrieve the actual XICS from the spectra
        retrieve_xics_from_raw_spectra(raw_files.combine(xics_to_extract))
        
        // IFF Present, load the trafoXMLs
        non_empty_trafoXMLs = trafoXMLs.ifEmpty(file("NO_TRAFO.trafoxml")).map({ it -> tuple("$it.baseName", it) })

        xics_and_trafos = retrieve_xics_from_raw_spectra.out
            .groupTuple(by: 0)
            .join(non_empty_trafoXMLs, remainder: true)
            .filter{ it[1] != null }
            .transpose(by: [1])
            .map({ it -> tuple(it[0], it[1], it[2] == null ? file("NULL___" + new java.util.Random().nextLong() + ".trafoXML") : it[2] ) })

        // Visualize in Plotly Plots
        if (params.export_visualization) {
            visualize_xics_via_plotly(
                xics_and_trafos.map({ it -> it[1] }).collect(),
                xics_and_trafos.map({ it -> it[2] }).collect()
            )
        }

    emit:
        xics_and_trafos.map({ it -> it[1] }).collect()  // HDF5 files
        xics_and_trafos.map({ it -> it[2] }).collect()  // Optional TrafoXML (otherwise empty files, with "NULL___" prefixed)
}


// Actual retrieval of the XICs using TRFP (Wrapper)
process retrieve_xics_from_raw_spectra {
    container "luxii/xic-extractor:latest"
    publishDir "${params.outdir}/extracted_xics/", mode:'copy'

    cpus 1
    memory "12 GB"
    stageInMode "link"

    input:
    tuple path(raw_file), path(queries_json)

    output:
    tuple val("${raw_file.baseName}"), path("${raw_file.baseName}.hdf5")

    """
    extract_via_fisher.py -raw $raw_file -query_csv $queries_json -out_hdf5 ${raw_file.baseName}.hdf5

    # TODO do this for Bruker!
    """
}


process visualize_xics_via_plotly {
    container "luxii/xic-extractor:latest"
    publishDir "${params.outdir}/extracted_xics_visualized/", mode:'copy'

    input:
    path(hdf5_xics)
    path(trafoXMLs_optional)

    output:
    tuple file("*.html"), file("*.png")
    """
    CONCAT_HDF5=""
    for file in $hdf5_xics
    do
        CONCAT_HDF5+="\$file,"
    done
    CONCAT_HDF5=\$(echo \$CONCAT_HDF5 | rev | cut -c2- | rev)

    CONCAT_TRAFO=""
    for file in $trafoXMLs_optional
    do
        CONCAT_TRAFO+="\$file,"
    done
    CONCAT_TRAFO=\$(echo \$CONCAT_TRAFO | rev | cut -c2- | rev)

    visualize_xic.py -input_hdf5s \$CONCAT_HDF5 -input_trafos \$CONCAT_TRAFO -group_regex '$params.group_regex' -outdir .
    """
}
