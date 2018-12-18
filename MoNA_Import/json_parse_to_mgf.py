#!/usr/bin/python


import sys
import getopt
import os
import requests
import json
import ming_fileio_library
import ming_spectrum_library
import cPickle as pickle

acceptable_instruments = ["DIRECT INFUSION NANOESI-ION TRAP", "LC-ESI-IT", "LC-ESI-ITFT", "LC-ESI-ITTOF", "LC-ESI-QFT", "LC-ESI-QIT", "LC-ESI-QQ", "LC-ESI-QTOF", "LC-ESI-TOF"]

def usage():
    print "<input json folder> <output_mgf> <output annotation table> <old library file>"


def parse_spectrum(spectrum):
    is_gnps = False
    is_hmdb = False
    is_insilico = False
    #Detecting if GNPS
    for tag in spectrum["tags"]:
        if tag["text"] == "gnps":
            is_gnps = True
            break
        if tag["text"] == "hmdb":
            is_hmdb = True
            break
        if tag["text"].find("virtual") != -1:
            is_insilico = True
            break
    if is_gnps:
        return None
    if is_hmdb:
        return None
    if is_insilico:
        return None


    inchi = spectrum["chemicalCompound"]["inchi"]
    compound_name = spectrum["chemicalCompound"]["names"][0]["name"]
    compound_id = spectrum["id"]
    adduct = ""
    instrument_type = ""
    precursor_mz = 0.0
    ionmode = ""
    cas_number = "N/A"

    for metadatum in spectrum["metaData"]:
        #print metadatum
        #print metadatum["name"]
        if metadatum["name"] == "precursortype" or metadatum["name"] == "precursor type":
            adduct = metadatum["value"]
        if metadatum["name"] == "instrument type":
            instrument_type = metadatum["value"]
        if metadatum["name"] == "precursor m/z":
            precursor_mz = float(metadatum["value"])
        if metadatum["name"] == "CAS":
            cas_number = metadatum["value"]

    #Validating Stuff
    if len(instrument_type) < 1:
        print "WITHOUTINST\t" + str(compound_id)
        return None

    if not instrument_type in acceptable_instruments:
        return None

    if len(adduct) < 1:
        print "WITHOUTADDUCT\t" + str(compound_id)
        return None

    if precursor_mz == 0.0:
        print "WITHOUTMZ\t" + str(compound_id)
        return None


    parsed_peaks = get_spectrum_peaks(spectrum)

    spectrum_obj = ming_spectrum_library.Spectrum("none", "1", "1", parsed_peaks, precursor_mz, 1)
    library_spectrum_obj = ming_spectrum_library.LibrarySpectrum(spectrum_obj)
    library_spectrum_obj.compound_name = "MoNA:" + str(compound_id) + " " + compound_name
    library_spectrum_obj.instrument = instrument_type
    library_spectrum_obj.ionsource = instrument_type
    library_spectrum_obj.inchi = inchi
    library_spectrum_obj.ionmode = ionmode
    library_spectrum_obj.acquisition = "isolated"
    library_spectrum_obj.collector = "MoNA"
    library_spectrum_obj.libraryquality = "3"
    library_spectrum_obj.CAS = cas_number
    library_spectrum_obj.pi = "MoNA"
    library_spectrum_obj.adduct = adduct
    return library_spectrum_obj


def get_spectrum_peaks(spectrum):
    peaks = spectrum["spectrum"].split(" ")
    parsed_peaks = []
    for peak in peaks:
        mass = float(peak.split(":")[0])
        intensity = float(peak.split(":")[1])
        parsed_peaks.append([mass, intensity])

    return parsed_peaks


def main():
    json_folder = sys.argv[1]
    output_mgf_filename = sys.argv[2]
    output_table_filename = sys.argv[3]

    input_old_library_filename = sys.argv[4]

    output_mgf_file = open(output_mgf_filename, "w")
    output_table_file = open(output_table_filename, "w")

    all_json_files = ming_fileio_library.list_files_in_dir(json_folder)

    #Current Library
    new_library = []

    for json_file in all_json_files:
        print json_file
        json_obj = json.loads(open(json_file, "r").read())
        for spectrum in json_obj:
            library_spectrum_obj = parse_spectrum(spectrum)

            if library_spectrum_obj != None:
                new_library.append(library_spectrum_obj)

    #Checking new library against the old one
    #old_library_list = ming_spectrum_library.load_mgf_file(input_old_library_filename)
    old_library_list = pickle.load(open(input_old_library_filename, 'rb'))

    output_spectrum_count = 0

    #Writing output table header
    output_table_file.write(ming_spectrum_library.LibrarySpectrum.get_gnps_library_creation_header() + "\n")

    for new_library_spectrum in new_library:
        max_similarity = 0.0
        for old_library_spectrum in old_library_list:
            old_mz = old_library_spectrum.mz
            new_lib_mz = new_library_spectrum.spectrum.mz

            mz_diff = abs(old_mz - new_lib_mz)
            if mz_diff > 0.1:
                continue

            cos_score = old_library_spectrum.cosine_spectrum(new_library_spectrum.spectrum, 1.0)
            max_similarity = max(max_similarity, cos_score)

        if max_similarity < 0.9:
            print new_library_spectrum.compound_name + "\t" + str(max_similarity)
            output_spectrum_count += 1
            new_library_spectrum.spectrum.scan = output_spectrum_count
            output_mgf_file.write(new_library_spectrum.get_mgf_string().encode('utf-8'))
            output_table_file.write(new_library_spectrum.get_gnps_library_creation_tsv_string(output_mgf_filename).encode('utf-8'))






if __name__ == "__main__":
    main()
