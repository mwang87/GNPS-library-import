import pandas as pd
import argparse
import json

import sys
sys.path.insert(0, 'ming_library')
import ming_fileio_library
import ming_spectrum_library
import ming_proteosafe_library
import ming_gnps_library


def parse_spectrum(spectrum_dict):
    peaks_json = json.loads(spectrum_dict["spectrum"])
    parsed_peaks = list(zip(peaks_json[0], peaks_json[1]))

    compound_id = spectrum_dict["Unnamed: 0"]
    
    precursor_mz = spectrum_dict["precursor_mz"]
    compound_name = spectrum_dict["name"]
    instrument_type = "Orbitrap"
    adduct = spectrum_dict["adduct"]
    smiles = spectrum_dict["smiles"]
    inchi = spectrum_dict["inchi"]
    smiles = spectrum_dict["smiles"]
    ionmode = spectrum_dict["polarity"]

    if ionmode == "positive":
        ionmode = "Positive"
    if ionmode == "negative":
        ionmode = "Negative"

    adduct = adduct[1:-2]

    spectrum_obj = ming_spectrum_library.Spectrum("none", "1", "1", parsed_peaks, precursor_mz, charge=1, ms_level=2)
    library_spectrum_obj = ming_spectrum_library.LibrarySpectrum(spectrum_obj)
    library_spectrum_obj.compound_name = compound_name.rstrip().replace("\n", "")
    library_spectrum_obj.instrument = instrument_type.rstrip().replace("\n", "")
    library_spectrum_obj.ionsource = "LC-ESI"
    library_spectrum_obj.inchi = inchi.rstrip().replace("\n", "")
    library_spectrum_obj.smiles = smiles
    library_spectrum_obj.ionmode = ionmode.rstrip().replace("\n", "")
    library_spectrum_obj.acquisition = "commercial"
    library_spectrum_obj.collector = "JGI:" + str(compound_id)
    library_spectrum_obj.libraryquality = "3"
    library_spectrum_obj.CAS = "N/A"
    library_spectrum_obj.pi = "Trent Northen"
    library_spectrum_obj.adduct = adduct.rstrip().replace("\n", "")

    return library_spectrum_obj

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("input_filename")
parser.add_argument("output_tsv")
parser.add_argument("output_mgf")

args = parser.parse_args()

print(args)

df = pd.read_csv(args.input_filename)

print(df.head())

input_list = df.to_dict(orient="records")

output_mgf_file = open(args.output_mgf, "w")
output_table_file = open(args.output_tsv, "w")

new_library_filtered = []
for input_dict in input_list:
    new_spectrum = parse_spectrum(input_dict)
    new_library_filtered.append(new_spectrum)

output_table_file.write(ming_spectrum_library.LibrarySpectrum.get_gnps_library_creation_header() + "\n")
output_spectrum_count = 0
for new_library_spectrum in new_library_filtered:
    output_spectrum_count += 1
    new_library_spectrum.spectrum.scan = output_spectrum_count
    output_mgf_file.write(new_library_spectrum.get_mgf_string())
    output_table_file.write(new_library_spectrum.get_gnps_library_creation_tsv_string(args.output_mgf).rstrip() + "\n")

output_table_file.close()
output_mgf_file.close()