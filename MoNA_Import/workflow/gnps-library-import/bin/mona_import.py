#!/usr/bin/python


import sys
import getopt
import os
import requests
import json
import zipfile
import shutil
import ming_fileio_library
import ming_spectrum_library
import ming_proteosafe_library
import ming_gnps_library

acceptable_instruments = ["DIRECT INFUSION NANOESI-ION TRAP", "LC-ESI-IT", "LC-ESI-ITFT", "LC-ESI-ITTOF", "LC-ESI-QFT", "LC-ESI-QIT", "LC-ESI-QQ", "LC-ESI-QTOF", "LC-ESI-TOF", "Linear Ion Trap", "ESI-QTOF", "ESI-QFT"]

def usage():
    print("<workflow xml> <input json folder> <output_mgf> <output annotation table> <output_script>")


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

    inchi = spectrum["compound"][0]["inchi"]
    compound_name = spectrum["compound"][0]["names"][0]["name"]
    compound_id = spectrum["id"]
    adduct = ""
    instrument_type = ""
    precursor_mz = 0.0
    ionmode = ""
    cas_number = "N/A"

    for metadatum in spectrum["metaData"]:
        if metadatum["name"] == "precursortype" or metadatum["name"] == "precursor type":
            adduct = metadatum["value"]
        if metadatum["name"] == "instrument type":
            instrument_type = metadatum["value"]
        if metadatum["name"] == "precursor m/z":
            precursor_mz = float(metadatum["value"])
        if metadatum["name"] == "CAS":
            cas_number = metadatum["value"]
        if metadatum["name"] == "ionization mode":
            ionmode = metadatum["value"]

    #Validating Stuff
    if len(instrument_type) < 1:
        print("WITHOUTINST\t" + str(compound_id))
        return None

    if not instrument_type in acceptable_instruments:
        print("BADINST\t" + str(instrument_type))
        return None

    if len(adduct) < 1:
        print("WITHOUTADDUCT\t" + str(compound_id))
        return None

    if precursor_mz == 0.0:
        print("WITHOUTMZ\t" + str(compound_id))
        return None


    parsed_peaks = get_spectrum_peaks(spectrum)

    spectrum_obj = ming_spectrum_library.Spectrum("none", "1", "1", parsed_peaks, precursor_mz, charge=1, ms_level=2)
    library_spectrum_obj = ming_spectrum_library.LibrarySpectrum(spectrum_obj)
    library_spectrum_obj.compound_name = compound_name.rstrip().replace("\n", "")
    library_spectrum_obj.instrument = instrument_type.rstrip().replace("\n", "")
    library_spectrum_obj.ionsource = "N/A"
    library_spectrum_obj.inchi = inchi.rstrip().replace("\n", "")
    library_spectrum_obj.ionmode = ionmode.rstrip().replace("\n", "")
    library_spectrum_obj.acquisition = "isolated"
    library_spectrum_obj.collector = "MoNA:" + str(compound_id)
    library_spectrum_obj.libraryquality = "3"
    library_spectrum_obj.CAS = cas_number.rstrip().replace("\n", "")
    library_spectrum_obj.pi = "MoNA"
    library_spectrum_obj.adduct = adduct.rstrip().replace("\n", "")

    library_spectrum_obj.MONAID = "MoNA:" + str(compound_id)

    return library_spectrum_obj


def get_spectrum_peaks(spectrum):
    peaks = spectrum["spectrum"].split(" ")
    parsed_peaks = []
    for peak in peaks:
        mass = float(peak.split(":")[0])
        intensity = float(peak.split(":")[1])
        parsed_peaks.append([mass, intensity])

    return parsed_peaks

def pulldown_mona_json(json_url, output_json_filename):
    cmd = "wget %s -O json.zip" % (json_url)
    os.system(cmd)

    with zipfile.ZipFile('json.zip', 'r') as myzip:
        if len(myzip.namelist()) == 1:
            myzip.extract(myzip.namelist()[0])
            shutil.move(os.path.basename(myzip.namelist()[0]), output_json_filename)
        else:
            exit(1)


def main():
    param_xml_filepath = sys.argv[1]
    output_folder = sys.argv[2]

    output_mgf_file = open(os.path.join(output_folder, "library_mgf.mgf"), "w")
    output_table_file = open(os.path.join(output_folder, "annotation_table.tsv"), "w")

    ###Reading XML file
    params_obj = ming_proteosafe_library.parse_xml_file(open(param_xml_filepath))
    json_library_filename = "library.json"
    pulldown_mona_json(params_obj["MONA_JSON_URL"][0], json_library_filename)

    new_library = []

    for spectrum in json.loads(open(json_library_filename, "r").read()):
        library_spectrum_obj = parse_spectrum(spectrum)

        if library_spectrum_obj != None:
            new_library.append(library_spectrum_obj)

    """Checking Old Library"""
    existing_gnps_library = ming_gnps_library.pulldown_library(params_obj["TARGET_LIBRARY_NAME"][0])

    new_library_filtered = []

    for new_library_spectrum in new_library:
        does_exist = False
        #print(new_library_spectrum.compound_name, new_library_spectrum.MONAID)

        for old_library_spectrum in existing_gnps_library:
            if old_library_spectrum["Compound_Name"].find(new_library_spectrum.MONAID) != -1:
                does_exist = True
                break
            if old_library_spectrum["Data_Collector"].find(new_library_spectrum.MONAID) != -1:
                does_exist = True
                break

        if not does_exist:
            new_library_filtered.append(new_library_spectrum)
        else:
            print("EXISTS")



    output_table_file.write(ming_spectrum_library.LibrarySpectrum.get_gnps_library_creation_header() + "\n")
    output_spectrum_count = 0
    for new_library_spectrum in new_library_filtered:
        output_spectrum_count += 1
        new_library_spectrum.spectrum.scan = output_spectrum_count
        output_mgf_file.write(new_library_spectrum.get_mgf_string())
        output_table_file.write(new_library_spectrum.get_gnps_library_creation_tsv_string("library_mgf.mgf"))




    #
    #
    #
    #
    # #Current Library
    # new_library = []
    #
    # for json_file in all_json_files:
    #     print json_file
    #     json_obj = json.loads(open(json_file, "r").read())
    #     for spectrum in json_obj:
    #         library_spectrum_obj = parse_spectrum(spectrum)
    #
    #         if library_spectrum_obj != None:
    #             new_library.append(library_spectrum_obj)
    #
    # #Checking new library against the old one
    # #old_library_list = ming_spectrum_library.load_mgf_file(input_old_library_filename)
    # old_library_list = pickle.load(open(input_old_library_filename, 'rb'))
    #
    # output_spectrum_count = 0
    #
    # #Writing output table header
    # output_table_file.write(ming_spectrum_library.LibrarySpectrum.get_gnps_library_creation_header() + "\n")
    #
    # for new_library_spectrum in new_library:
    #     max_similarity = 0.0
    #     for old_library_spectrum in old_library_list:
    #         old_mz = old_library_spectrum.mz
    #         new_lib_mz = new_library_spectrum.spectrum.mz
    #
    #         mz_diff = abs(old_mz - new_lib_mz)
    #         if mz_diff > 0.1:
    #             continue
    #
    #         cos_score = old_library_spectrum.cosine_spectrum(new_library_spectrum.spectrum, 1.0)
    #         max_similarity = max(max_similarity, cos_score)
    #
    #     if max_similarity < 0.9:
    #         print new_library_spectrum.compound_name + "\t" + str(max_similarity)
    #         output_spectrum_count += 1
    #         new_library_spectrum.spectrum.scan = output_spectrum_count
    #         output_mgf_file.write(new_library_spectrum.get_mgf_string().encode('utf-8'))
    #         output_table_file.write(new_library_spectrum.get_gnps_library_creation_tsv_string(output_mgf_filename).encode('utf-8'))
    #
    #




if __name__ == "__main__":
    main()
