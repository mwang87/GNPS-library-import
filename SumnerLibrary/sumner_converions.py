#!/usr/bin/python


import sys
import getopt
import requests
import requests_cache
import json
import os

requests_cache.install_cache('demo_cache', allowable_codes=(200, 404))

def inchikey_to_inchi_smiles_pubchem(inchikey):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/%s/JSON" % (inchikey)
    #print(url)
    inchi = "N/A"
    smiles = "N/A"

    try:
        r = requests.get(url)
        if r.status_code == 200:
            json_obj = json.loads(r.text)
            for compound_property in json_obj["PC_Compounds"][0]["props"]:
                if compound_property["urn"]["label"] == "InChI":
                    inchi = compound_property["value"]["sval"]
                if compound_property["urn"]["label"] == "SMILES" and compound_property["urn"]["name"] == "Canonical":
                    smiles = compound_property["value"]["sval"]

    except KeyboardInterrupt:
        raise
    except:
        raise
        return inchi, smiles
    return inchi, smiles

def inchikey_to_inchi_chemspider(inchikey):
    url = "http://www.chemspider.com/InChI.asmx/InChIKeyToInChI"
    payload = {'inchi_key': inchikey}
    inchi = "N/A"
    try:
        r = requests.get(url, params=payload)
        if r.status_code == 200:
            inchi = xmltodict.parse(r.text)["string"]["#text"]
    except KeyboardInterrupt:
        raise
    except:
        return inchi
    return inchi

def inchi_to_smiles_chemspider(inchi):
    url = "http://www.chemspider.com/InChI.asmx/InChIToSMILES"
    payload = {'inchi': inchi}
    smiles = "N/A"
    try:
        r = requests.get(url, params=payload)
        if r.status_code == 200:
            smiles = xmltodict.parse(r.text)["string"]["#text"]
    except KeyboardInterrupt:
        raise
    except:
        return smiles
    return smiles

def main():
    txt_file = open(sys.argv[1], "r")
    mgf_filename = sys.argv[2]
    mgf_file = open(sys.argv[2], "w")
    batch_file = open(sys.argv[3], "w")


    acceptable_ionization = set(["ESI", "APCI"])
    acceptable_instruments = set([])

    acceptable_instruments = set(["ESI-IT-MS/MS", "ESI-QqIT-MS/MS", "ESI-QqQ-MS/MS", "ESI-QqTOF-MS/MS", "FAB-EBEB", "LC-APPI-QQ", "LC-ESI-IT", "LC-ESI-ITFT", "LC-ESI-ITTOF", "LC-ESI-QIT", "LC-ESI-QQ", "LC-ESI-QTOF"  ])

    peptide = "*..*"
    smiles = "N/A"
    inchi = "N/A"
    pepmass = "0.0"
    title = ""
    instrument = ""
    compound_name = ""
    peaks = []
    retentiontime = ""
    peaks_start = 0;
    exactmass = "0"
    cas_number = "N/A"
    adduct = "[M+H]"
    spectrum_level = 0
    ionization_mode = ""
    collision_energy = ""

    read_peaks = False

    scan_number = 1

    #Writing Batch Headers
    batch_file.write("FILENAME\tSEQ\tCOMPOUND_NAME\tMOLECULEMASS\tINSTRUMENT\tIONSOURCE\tEXTRACTSCAN\t")
    batch_file.write("SMILES\tINCHI\tINCHIAUX\tCHARGE\tIONMODE\tPUBMED\tACQUISITION\tEXACTMASS\tDATACOLLECTOR\t")
    batch_file.write("ADDUCT\tINTEREST\tLIBQUALITY\tGENUS\tSPECIES\tSTRAIN\tCASNUMBER\tPI\n")

    for line in txt_file:
        if line[:5] == "Name:":
            compound_name = line.strip()[len("Name: "):].strip()
            #print line.rstrip()

        if line.find("InstName") != -1:
            instrument = line.rstrip().split(":")[1]

        if line.find("MSMS") != -1:
            ms_level = line.split(":")[1].strip()

        if line.find("IonPolarity") != -1:
            ionization_mode = line.split(":")[1].strip()
            if ionization_mode == "neg":
                adduct = "M-H"
                ionization_mode = "Negative"
            else :
                adduct = "M+H"
                ionization_mode = "Positive"

        if line.find("Smiles") != -1:
            smiles = line.split(":")[1].strip()

        if line[:len("InChI:")] == "InChI:":
            inchi = line.split(":")[1].strip()

        if line.find("PreIon: ") != -1:
            pepmass = line[len("PreIon: "):].rstrip()

        if line.find("CAS:") != -1:
            cas_number = line[len("CAS: "):].rstrip()

        if line.find("ColEnergy:") != -1:
            collision_energy = line[len("ColEnergy: "):].rstrip()

        if line.find("Num Peaks: ") != -1:
            peaks = []
            read_peaks = True
            continue

        if len(line.strip()) < 1:
            if pepmass == "0.0":
                cas_number = "N/A"
                smiles = "N/A"
                inchi = "N/A"
                pepmass = "0.0"
                read_peaks = False
                continue
            #print(compound_name, pepmass)

            #End of spectrum, writing spectrum
            spectrum_string = ""
            spectrum_string += "BEGIN IONS\n"
            spectrum_string += "SEQ=" + peptide + "\n"
            spectrum_string += "PEPMASS=" + pepmass + "\n"
            spectrum_string += "SMILES=" + smiles + "\n"
            spectrum_string += "INCHI=" + inchi + "\n"
            spectrum_string += "SOURCE_INSTRUMENT=" + instrument + "\n"
            spectrum_string += "NAME=" + compound_name + "\n"
            spectrum_string += "ORGANISM=SUMNER\n"
            spectrum_string += "SCANS=" + str(scan_number) + "\n"

            for peak in peaks:
                spectrum_string += peak + "\n"


            peaks = []
            spectrum_string += "END IONS\n"
            #print spectrum_string
            mgf_file.write(spectrum_string)
            read_peaks = False

            #writing batch file
            batch_file.write(mgf_filename + "\t")
            batch_file.write(peptide + "\t")
            batch_file.write(compound_name + " - " + collision_energy + "eV\t")
            batch_file.write(pepmass + "\t")
            batch_file.write(instrument + "\t")
            batch_file.write("LC-ESI" + "\t")
            batch_file.write(str(scan_number) + "\t")
            batch_file.write(smiles + "\t")
            batch_file.write(inchi + "\t")
            batch_file.write("N/A" + "\t")
            batch_file.write("1" + "\t")
            batch_file.write(ionization_mode + "\t")
            batch_file.write("N/A" + "\t")
            batch_file.write("Isolated" + "\t")
            batch_file.write("0" + "\t")
            batch_file.write("Sumner" + "\t")
            batch_file.write(adduct + "\t")
            batch_file.write("N/A" + "\t")
            batch_file.write("3" + "\t")
            batch_file.write("N/A" + "\t")
            batch_file.write("N/A" + "\t")
            batch_file.write("N/A" + "\t")
            batch_file.write(cas_number + "\t")
            batch_file.write("Sumner" + "\n")


            scan_number += 1

            cas_number = "N/A"
            smiles = "N/A"
            inchi = "N/A"
            pepmass = "0.0"
            collision_energy = ""

        if read_peaks == True:
            peaks_splits = line.rstrip().split()
            peaks_joined = zip(peaks_splits[0::2], peaks_splits[1::2])
            for peak_join in peaks_joined:
                peaks.append(" ".join(peak_join))
            #peaks += peaks_joined

    return 0;




if __name__ == "__main__":
    main()
