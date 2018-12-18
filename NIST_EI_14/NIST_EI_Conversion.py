#!/usr/bin/python


import sys
import getopt
import requests
import requests_cache
import xmltodict
import json

requests_cache.install_cache('demo_cache', allowable_codes=(200, 404))
#requests_cache.install_cache('redis', backend='redis', allowable_codes=(200, 404))


def usage():
    print "<input txt> <output mgf> <output batch file> <inchikey mapping to inchi file> <inchikey mapping to smiles>"
    print "Takes NIST MSP file to convert to MGF and batch file"


def load_inchikey_mapping(inchi_mapping, smiles_mapping):
    inchikey_to_inchi_mapping = {}
    inchikey_to_smiles_mapping = {}

    try:
        for line in open(inchi_mapping):
            inchi_key = line.split("\t")[0]
            inchi_string = line.split("\t")[1].rstrip()
            inchikey_to_inchi_mapping[inchi_key] = inchi_string
    except:
        inchikey_to_inchi_mapping = {}

    try:
        for line in open(smiles_mapping):
            inchi_key = line.split("\t")[0]
            smiles_string = line.split("\t")[1].rstrip()
            inchikey_to_smiles_mapping[inchi_key] = smiles_string
    except:
        inchikey_to_smiles_mapping = {}

    return inchikey_to_inchi_mapping, inchikey_to_smiles_mapping

def inchikey_to_inchi_smiles_pubchem(inchikey):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/%s/JSON" % (inchikey)
    print(url)
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

def inchikey_to_inchi(inchikey):
    url = "https://cactus.nci.nih.gov/chemical/structure/%s/stdinchi" % (inchikey)
    print(url)
    inchi = "N/A"
    try:
        r = requests.get(url)
        print("STATUS", r.status_code)
        if r.status_code == 200:
            inchi = r.text
    except KeyboardInterrupt:
        raise
    except:
        return inchi
    return inchi


def inchikey_to_smiles(inchikey):
    url = "https://cactus.nci.nih.gov/chemical/structure/%s/smiles" % (inchikey)
    smiles = "N/A"
    try:
        r = requests.get(url)
        if r.status_code == 200:
            smiles = r.text
    except KeyboardInterrupt:
        raise
    except:
        return smiles
    return smiles

def main():
    usage()

    txt_file = open(sys.argv[1], "r")
    mgf_filename = sys.argv[2]
    mgf_file = open(sys.argv[2], "w")
    batch_file = open(sys.argv[3], "w")

    inchikey_to_inchi_mapping = {}
    inchikey_to_smiles_mapping = {}
    try:
        inchikey_to_inchi_mapping, inchikey_to_smiles_mapping = load_inchikey_mapping(sys.argv[4], sys.argv[5])
    except:
        print("No mapping")



    print sys.argv[1]

    peptide = "*..*"
    smiles = "N/A"
    inchi = "N/A"
    pepmass = ""
    title = ""
    instrument = "EI"
    compound_name = ""
    peaks = []
    retentiontime = ""
    ion_mode = "EI"
    peaks_start = 0;
    exactmass = "0"
    cas_number = "N/A"
    adduct = "EI"
    spectrum_level = 0
    ionization_mode = "EI"
    nist_no = " "

    read_peaks = False

    scan_number = 1

    #Writing Batch Headers
    batch_file.write("FILENAME\tSEQ\tCOMPOUND_NAME\tMOLECULEMASS\tINSTRUMENT\tIONSOURCE\tEXTRACTSCAN\t")
    batch_file.write("SMILES\tINCHI\tINCHIAUX\tCHARGE\tIONMODE\tPUBMED\tACQUISITION\tEXACTMASS\tDATACOLLECTOR\t")
    batch_file.write("ADDUCT\tINTEREST\tLIBQUALITY\tGENUS\tSPECIES\tSTRAIN\tCASNUMBER\tPI\n")

    for line in txt_file:

        if line.find("Name:") != -1:
            compound_name = line.strip()[len("Name: "):]

        if line.find("NIST#:") != -1:
            nist_no = line[len("NIST#: "):].rstrip()

        if line.find("ExactMass: ") != -1:
            pepmass = line[len("ExactMass: "):].rstrip()
            exactmass = line[len("ExactMass: "):].rstrip()

        if line.find("CAS#: ") == 0:
            cas_number = line.split(";")[0][len("CAS#: "):].rstrip()
            nist_no = line.split(";")[1].strip().replace("NIST#: ", "")

        if line.find("Num Peaks: ") != -1:
            peaks = []
            read_peaks = True
            continue

        if line.find("InChIKey: ") != -1:
            inchi_key = line[len("InChIKey: "):].rstrip()
            #print("INCHI\t" + inchi_key)
            if inchi_key in inchikey_to_inchi_mapping:
                inchi = inchikey_to_inchi_mapping[inchi_key]
            if inchi_key in inchikey_to_smiles_mapping:
                smiles = inchikey_to_smiles_mapping[inchi_key]
            #if len(inchi) < 4:
            #    inchi, smiles = inchikey_to_inchi_smiles_pubchem(inchi_key)
            #    print(inchi, smiles)

            # inchi = inchikey_to_inchi_chemspider(inchi_key)
            # smiles = inchi_to_smiles_chemspider(inchi)
            # if inchi == "N/A":
            #     print("CACTUS LOOKUP", inchi_key)
            #     inchi = inchikey_to_inchi(inchi_key)
            #     smiles = inchikey_to_smiles(inchi_key)
            # smiles = smiles.replace("\n","")
            # inchi = inchi.replace("\n","")



        if len(line.strip()) < 1:
            #End of spectrum, writing spectrum
            spectrum_string = ""
            spectrum_string += "BEGIN IONS\n"
            spectrum_string += "SEQ=" + peptide.rstrip() + "\n"
            spectrum_string += "PEPMASS=" + pepmass.rstrip() + "\n"
            spectrum_string += "SMILES=" + smiles.rstrip() + "\n"
            spectrum_string += "INCHI=" + inchi.rstrip() + "\n"
            spectrum_string += "SOURCE_INSTRUMENT=" + instrument.rstrip() + "\n"
            spectrum_string += "NAME=" + "NIST:" + nist_no.rstrip() + " " + compound_name.rstrip() + "\n"
            spectrum_string += "ORGANISM=NIST\n"
            spectrum_string += "SCANS=" + str(scan_number).rstrip() + "\n"

            for peak in peaks:
                peaks_cleaned = peak.strip()
                if len(peaks_cleaned) > 2:
                    spectrum_string += peaks_cleaned + "\n"


            peaks = []
            spectrum_string += "END IONS\n"
            #print spectrum_string
            mgf_file.write(spectrum_string)
            read_peaks = False

            #writing batch file
            batch_file.write(mgf_filename.rstrip() + "\t")
            batch_file.write(peptide.rstrip() + "\t")
            batch_file.write(compound_name.rstrip() + "\t")
            batch_file.write(pepmass.rstrip() + "\t")
            batch_file.write(instrument.rstrip() + "\t")
            batch_file.write(ionization_mode.rstrip() + "\t")
            batch_file.write(str(scan_number).rstrip() + "\t")
            batch_file.write(smiles.rstrip() + "\t")
            batch_file.write(inchi.rstrip() + "\t")
            batch_file.write("N/A" + "\t")
            batch_file.write("1" + "\t")
            batch_file.write(ion_mode.rstrip() + "\t")
            batch_file.write("N/A" + "\t")
            batch_file.write("Isolated" + "\t")
            batch_file.write("0" + "\t")
            batch_file.write("NIST" + "\t")
            batch_file.write(adduct.rstrip() + "\t")
            batch_file.write("N/A" + "\t")
            batch_file.write("3" + "\t")
            batch_file.write("N/A" + "\t")
            batch_file.write("N/A" + "\t")
            batch_file.write("N/A" + "\t")
            batch_file.write(cas_number.rstrip() + "\t")
            batch_file.write("NIST" + "\n")


            scan_number += 1

            cas_number = "N/A"
            smiles = "N/A"
            inchi = "N/A"

            print(scan_number)
        if read_peaks == True:
            peaks += line.rstrip().split(";")
            #peaks.append(line.rstrip())

    return 0;




if __name__ == "__main__":
    main()
