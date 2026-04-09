# mzkitcpp
mass spectrometry utility functions and processing algorithms written in C++, exposed in R package. 

This package has a direct connection to [maven_core](https://github.com/eugenemel/maven_core) C++ library, and is the preferred library for wrapping [maven_core](https://github.com/eugenemel/maven_core) functionality into R. [maven_core](https://github.com/eugenemel/maven_core) implements all non-gui functionality associated with the [MAVEN GUI](https://github.com/eugenemel/maven/releases/latest) (see [maven repo](https://github.com/eugenemel/maven) for source code)

# Installation
Execute the following command in an R console:
```
remotes::install_github("calico/mzkitcpp", force=TRUE, build_vignettes=TRUE)
```
# Version Changelog

- **1.0.0** initial public release
- **1.1.0** `rsecprofiler`-related functionality
- **1.2.0** hrms qc app functionality

# Functions
- **DI_encoded_search_params**: Encodes DI search parameters into a formatted string.
- **DI_file_info**: Extracts file information and scan statistics from DI mzML files.
- **DI_ms1_and_ms2_intensity**: Extracts MS1 precursor and MS2 fragment intensities for specified m/z values from DI samples.
- **DI_ms1_range_intensity**: Extracts summed MS1 intensities within specified precursor m/z ranges.
- **DI_ms2_range_id**: Identifies which MS2 isolation window range each precursor m/z belongs to.
- **DI_ms2_ranges**: Extracts MS2 isolation window ranges from a DI mzML file.
- **DI_ms3_intensity**: Extracts MS3 fragment intensities for specified MS1 precursor, MS2 precursor, and MS3 fragment m/z values (typically used for triacylglyceride analysis on ID-X tribrid instruments).
- **DI_ms3_peakgroups_and_peaks**: Converts MS3 search results to peakgroup and peak table format for database storage.
- **DI_ms3_targets**: Extracts MS3 target information (MS1 and MS2 precursor m/z pairs) from an mzML file.
- **DI_ms_intensity**: Extracts intensities for specified m/z values at a given MS level from DI samples.
- **DI_peakgroups_and_peaks**: Converts DI search results to peakgroup and peak table format for database storage.
- **DI_pipeline**: Main pipeline for DI MS2 compound identification and quantification.
- **DI_pipeline_ms3_search**: Pipeline for MS3 analysis, specifying both MS1 and MS2 precursor m/z values (typically for triacylglycerides on ID-X tribrid instruments).
- **DI_search_lib**: Searches a spectral library against DI MS2 data for compound identification.
- **DI_slice_library**: Partitions a spectral library by MS2 isolation window ranges for efficient searching.
- **DI_summarized_compounds**: Summarizes compound identifications using parsimony and fragment evidence.
- **DI_unslice_library**: Reverses library slicing to restore original library structure.
- **ISO_isotope_matrices**: Generates isotope distribution matrices for specified compounds.
- **adductize_exact_mass**: Calculates precursor m/z from an exact mass and adduct.
- **adductize_formula**: Calculates precursor m/z from a molecular formula and adduct.
- **adductize_peptide**: Calculates precursor m/z from a peptide sequence and adduct.
- **envelope_dist_peptide**: Calculates the isotopic envelope distribution for a peptide sequence.
- **exact_mass**: Calculates the exact (monoisotopic) mass from a molecular formula.
- **exact_mass_peptide**: Calculates the exact (monoisotopic) mass from a peptide sequence.
- **expand_redundant_fragments_ms2_lib**: Expands MS2 library entries containing redundant fragments into separate rows.
- **expand_redundant_fragments_ms3_lib**: Expands MS3 library entries containing redundant fragments into separate rows.
- **export_msp_lipids_library**: Writes a spectral library DataFrame to an MSP format file.
- **findFragPairsGreedyMz**: (Deprecated) Matches library and experimental fragment m/z values using greedy algorithm.
- **find_duplicate_peaks**: Identifies duplicate peaks between original and re-picked peak tables.
- **get_background_subtracted_scan_data**: Retrieves scan data with background subtraction applied.
- **get_consensus_spectrum**: Creates a consensus MS2 spectrum from multiple scans within a sample.
- **get_eic**: Extracts an extracted ion chromatogram (EIC) for a specified m/z and RT range.
- **get_last_rts**: Retrieves the final retention time from each sample file.
- **get_ms2_scans**: Extracts MS2 scan information from sample files.
- **get_multi_file_consensus_spectrum**: Creates a consensus MS2 spectrum from scans across multiple sample files.
- **get_scan_data**: Retrieves m/z and intensity data for specified scan numbers.
- **get_scan_metadata**: Retrieves metadata (RT, precursor m/z, MS level, etc.) for specified scans.
- **get_standard_rts**: Retrieves retention times for internal standards across sample files.
- **groups_to_msp_library**: Converts peakgroup/peak tables to MSP spectral library format.
- **import_msp_lipids_library**: Imports an MSP format spectral library into a DataFrame.
- **lipidmaps_2020_compound_names**: Processes compound names using LipidMaps 2020 nomenclature.
- **loop_injections_to_msp_library**: Converts loop injection data to MSP spectral library format.
- **maldesi_create_modified_mzML**: Creates a modified mzML file with adjusted scan parameters for MALDESI data.
- **maldesi_isotopic_envelope_finder**: Identifies isotopic envelopes in MALDESI imaging mass spectrometry data.
- **maldesi_search**: Searches for compounds in MALDESI imaging data using exact mass or adduct-based matching.
- **mark_fragments_ms2_lib**: Marks specific fragments in an MS2 spectral library (e.g., as diagnostic ions).
- **match_id_candidates_rcpp**: (Deprecated) Matches compound identification candidates based on m/z.
- **merge_split_groups**: Identifies and merges peakgroups that were incorrectly split during peak detection.
- **monoiosotopic_mass**: Calculates monoisotopic masses and adducted m/z values for molecular formulas.
- **mzk_get_isotope_parameters**: Encodes isotope search parameters into a formatted string.
- **mzk_get_lipid_parameters**: Encodes lipid search parameters into a formatted string.
- **name_summaries**: Parses and standardizes compound names into structured components.
- **peptide_sequence_to_formula**: Converts a peptide amino acid sequence to its molecular formula.
- **precursor_mass**: Calculates precursor m/z values from molecular formulas and adducts.
- **predict_formula**: Predicts possible molecular formulas from a measured m/z value.
- **qqq_peaks**: Extracts and quantifies peaks from triple quadrupole (QQQ) mass spectrometry data.
- **record_set_to_msp_library**: Converts formatted metabolite records to MSP spectral library format.
- **reextract_peaks**: Re-extracts peak intensities from raw data using updated parameters.
- **repick_peaks**: Re-picks peak boundaries and quantifies peaks using updated parameters.
- **smoothed_series**: Applies smoothing algorithms to time series data.
- **update_rts**: Updates retention times based on an alignment key.
