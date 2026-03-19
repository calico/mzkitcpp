# mzkitcpp
mass spectrometry utility functions and processing algorithms written in C++, exposed in R package. 
This package has a direct connection to `maven_core` C++ library, and is the preferred library for
wrapping `maven_core` functionality into R.

# Installation
Execute the following command in an R console:
```
remotes::install_github("calico/mzkitcpp", force=TRUE, build_vignettes=TRUE)
```

# Functions
- **DI_encoded_search_params**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_file_info**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_ms1_and_ms2_intensity**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_ms1_range_intensity**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_ms2_range_id**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_ms2_ranges**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_ms3_intensity**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_ms3_peakgroups_and_peaks**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_ms3_targets**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_ms_intensity**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_peakgroups_and_peaks**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_pipeline**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_pipeline_ms3_search**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_search_lib**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_slice_library**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_summarized_compounds**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **DI_unslice_library**: Direct infusion mass spectrometry (DI) analysis pipeline function for MS/MS-all type DIMS analysis.
- **ISO_isotope_matrices**: R-facing functions
- **adductize_exact_mass**: Given an exact mass and an adduct name, return the precursor m/z.
- **adductize_formula**: Given a molecular formula and an adduct name, return the precursor m/z.
- **adductize_peptide**: Given a peptide sequence and an adduct name, return the precursor m/z.
- **envelope_dist_peptide**: Given a peptide sequence, return the exact mass.
- **exact_mass**: Given a molecular formula, return the exact mass.
- **exact_mass_peptide**: Given a peptide sequence, return the exact mass.
- **expand_redundant_fragments_ms2_lib**: Given an MS2 library data frame containing redundant fragments, expand these redundant fragments into separate rows, with all other columns duplicated.
- **expand_redundant_fragments_ms3_lib**: Given an MS3 library data frame containing redundant fragments, expand' these redundant fragments into separate rows, with all other columns duplicated.
- **export_msp_lipids_library**: Given a DataFrame that originally derived from an imported msp library, write the contents of the DataFrame to an msp file.
- **findFragPairsGreedyMz**: for debugging static void printArr(vector<pair < float,pair < int, int > > > arr);  @param library_peaks: list of m/zs corresponding to fragment spectra, ordered in increasing order @param exp_peaks: list of m/zs correspond to experimental spectra, ordered in increasing m/z @param maxMzDiffInAmu: maximum m/z difference between two fragments allowable for matching @param min_intensity_frag_exp: smallest allowable intensity fraction of max intensity for an experimental fragment.
- **find_duplicate_peaks**: Given a merged dataframe of a peaks-formatted data frames of original peaks and re-picked peaks, return a list of peaks that are sufficiently similar to both.
- **get_background_subtracted_scan_data**: Assume input scan number counts from 1 instead of 0 Scan information
- **get_consensus_spectrum**: Mass spectrometry processing utility.
- **get_eic**: Mass spectrometry processing utility.
- **get_last_rts**: Mass spectrometry processing utility.
- **get_ms2_scans**: samples_table is the data from the mzrollDB samples table.
- **get_multi_file_consensus_spectrum**: Mass spectrometry processing utility.
- **get_scan_data**: Assume input scan number counts from 1 instead of 0 Scan information
- **get_scan_metadata**: Different scans may be different samples, so this must match exactly.
- **get_standard_rts**: Mass spectrometry processing utility.
- **groups_to_msp_library**: samples_table is the data from the mzrollDB samples table.
- **import_msp_lipids_library**: Given an msp library, return a DataFrame containing relevant summary information output debugging
- **lipidmaps_2020_compound_names**: Given an msp library, return a DataFrame containing relevant summary information
- **loop_injections_to_msp_library**: metabolite_data is a data table that contains a path to a loop injection file, along with other columns with key information write information to msp file Ensure that this compound will not be recomputed again Determine which loop files may be bad close msp file, free memory
- **maldesi_create_modified_mzML**: corresponding sections of the mzML file adjusted appropriately, so that mzML parsing programs can interpret the data correctly.
- **maldesi_isotopic_envelope_finder**: clean up
- **maldesi_search**: If a set of adducts are provided as an input parameter, this value is treated as an exact mass, from which theoretical m/z values are computed.
- **mark_fragments_ms2_lib**: Given an MS2 library data frame containing redundant fragments, expand these redundant fragments into separate rows, with all other columns duplicated.
- **match_id_candidates_rcpp**: This function is a reimplementation of the matching part of clamdb<matching.
- **merge_split_groups**: Given a dataframe sorted by (groupId, sampleId), identify peakgroups that were split, and should be merged together (based on m/z and RT tolerance).
- **monoiosotopic_mass**: Given a vector of string formulas, return a vector of compound monoisotopic masses.
- **mzk_get_isotope_parameters**: Given various parameter values, return a formatted string as produced by IsotopeParameters::encode() isotope enumeration isotope extraction only consider max natural abundance error if the isIgnoreNaturalAbundance flag is TRUE - from MAVEN 1 diff iso specific
- **mzk_get_lipid_parameters**: Given various source information, return a formatted string as produced by LCLipidSearchParameters::encode()
- **name_summaries**: Given a vector of string formulas, return a vector of compound monoisotopic masses.
- **peptide_sequence_to_formula**: Convert a peptide sequence to its molecular formula string.
- **precursor_mass**: Given a vector of string formulas, return a vector of compound monoisotopic masses.
- **predict_formula**: Given a vector of string formulas, return a vector of compound monoisotopic masses.
- **qqq_peaks**: rt types quant types smoothed quant types Parameters always need to compute bounds to pass back complete information Samples rt types quant types smoothed quant types rt types quant types smoothed quant types rt types quant types smoothed quant types rt types quant types smoothed quant types
- **record_set_to_msp_library**: metabolite_data is a data table that contains formatted metabolite records information (for use with DBv3).
- **reextract_peaks**: clean up
- **repick_peaks**: @param peaks_to_repick: DataFrame with instructions for repicking.
- **smoothed_series**: Given a data vector and various other parameters, perform smoothing.
- **update_rts**: Mass spectrometry processing utility.
