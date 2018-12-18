#!/usr/bin/python

###
# This file take as input a concatonation of the txt files
# Then spits out the mgf and batch file
###

import sys
import getopt
import ming_gnps_library

def get_list_of_massbank_accessions_already_existing():
    all_library_spectra = []
    all_library_spectra += ming_gnps_library.pulldown_library("MASSBANK")
    all_library_spectra += ming_gnps_library.pulldown_library("MASSBANKEU")

    set_of_accessions = set()
    for library_spec in all_library_spectra:
        #print library_spec["Compound_Name"]
        #print library_spec
        if library_spec["Compound_Name"].find("Massbank") == -1:
            continue
        accession = library_spec["Compound_Name"].split(" ")[0].split(":")[1]
        #print accession
        set_of_accessions.add(accession)

    return set_of_accessions

def usage():
    print "<input txt> <output mgf> <output batch file>"




def main():
    usage()

    txt_file = open(sys.argv[1], "r")
    mgf_filename = sys.argv[2]
    mgf_file = open(sys.argv[2], "w")
    batch_file = open(sys.argv[3], "w")

    batch_file.write("FILENAME\tSEQ\tCOMPOUND_NAME\tMOLECULEMASS\tINSTRUMENT\tIONSOURCE\tEXTRACTSCAN\t")
    batch_file.write("SMILES\tINCHI\tINCHIAUX\tCHARGE\tIONMODE\tPUBMED\tACQUISITION\tEXACTMASS\tDATACOLLECTOR\t")
    batch_file.write("ADDUCT\tINTEREST\tLIBQUALITY\tGENUS\tSPECIES\tSTRAIN\tCASNUMBER\tPI\n")


    print sys.argv[1]

    set_of_existing_accessions = get_list_of_massbank_accessions_already_existing()

    acceptable_instruments = set(["ESI-IT-MS/MS", "ESI-QqIT-MS/MS", "ESI-QqQ-MS/MS", "ESI-QqTOF-MS/MS", "FAB-EBEB", "LC-APPI-QQ", "LC-ESI-IT", "LC-ESI-ITFT", "LC-ESI-ITTOF", "LC-ESI-QIT", "LC-ESI-QQ", "LC-ESI-QTOF", "ESI-ITFT", "UPLC-ESI-QTOF", "LC-ESI-Q", "LC-ESI-QFT", "ESI-FTICR", "HPLC-ESI-TOF", "LC-APCI-Q", "APCI-ITFT"])

    peptide = "*..*"
    smiles = ""
    inchi = ""
    pepmass = ""
    title = ""
    instrument = ""
    compound_name = ""
    peaks = []
    retentiontime = ""
    ion_mode = ""
    peaks_start = 0;
    exactmass = "0"
    cas_number = ""
    adduct = "[M+H]"
    spectrum_level = 0
    ion_mode = "Positive"

    scan_number = 1

    for line in txt_file:
        #writing out spectrum
        if line.find("//") != -1:
            write_spectrum = 1

            #print instrument

            if title in set_of_existing_accessions:
                print "Accession Exists"
                peaks = []
                write_spectrum = 0

            if instrument in acceptable_instruments:
                x = 1
                #if ion_mode == "NEGATIVE":
                    #x = 1
                #    peaks = []
                #    write_spectrum = 0
                #else:
                    #continue
            else:
                print "Unacceptalbe Instrument " + title + "\t" + instrument
                peaks = []
                write_spectrum = 0

            if spectrum_level != 2:
                #print "WRONG LEVEL"
                write_spectrum = 0

            if write_spectrum == 1:
                compound_name = compound_name[:-1]

                if len(compound_name) > 200:
                    #print title
                    #print compound_name
                    compound_name = compound_name[:200]
                    #print compound_name


                #looking at adducts
                if adduct == "[M+H]+":
                    ion_mode = "Positive"
                elif adduct == "[M-H]-":
                    ion_mode = "Negative"
                else:
                    print("Cannot find adduct", adduct)

                if len(pepmass) == 0:
                    resolved = 0
                    if len(exactmass) > 1 and adduct == "[M+H]+":
                        pepmass =  str(float(exactmass) + 1.007825)
                        resolved = 1
                    if len(exactmass) > 1 and adduct == "M+":
                        pepmass =  exactmass
                        resolved = 1
                    if len(exactmass) > 1 and adduct == "[M+H-H2O]+":
                        pepmass =  str(float(exactmass) + 1.007825 - 18.010565)
                        resolved = 1
                    if len(exactmass) > 1 and adduct == "[M+H-(C12H20O9)]+":
                        pepmass =  str(float(exactmass) + 1.007825 - 308.110735)
                        resolved = 1

                    if resolved == 0:
                        print "FUCK THIS SHIT: " + title

                spectrum_string = ""
                spectrum_string += "BEGIN IONS\n"
                spectrum_string += "SEQ=" + peptide + "\n"
                spectrum_string += "PEPMASS=" + pepmass + "\n"
                spectrum_string += "SMILES=" + smiles + "\n"
                spectrum_string += "INCHI=" + inchi + "\n"
                spectrum_string += "SOURCE_INSTRUMENT=" + instrument + "\n"
                spectrum_string += "NAME=" + "Massbank:" + title + " " + compound_name + "\n"
                spectrum_string += "ORGANISM=Massbank\n"
                spectrum_string += "RTINSECONDS=" + retentiontime + "\n"
                spectrum_string += "SCANS=" + str(scan_number) + "\n"


                batch_file.write(mgf_filename + "\t")
                batch_file.write(peptide + "\t")
                batch_file.write("MassbankEU:" + title + " " + compound_name + "\t")
                batch_file.write(pepmass + "\t")
                batch_file.write(instrument + "\t")
                batch_file.write("ESI" + "\t")
                batch_file.write(str(scan_number) + "\t")
                batch_file.write(smiles + "\t")
                batch_file.write(inchi + "\t")
                batch_file.write("N/A" + "\t")
                batch_file.write("1" + "\t")
                batch_file.write(ion_mode + "\t")
                batch_file.write("N/A" + "\t")
                batch_file.write("Isolated" + "\t")
                batch_file.write("0" + "\t")
                batch_file.write("Massbank" + "\t")
                batch_file.write(adduct + "\t")
                batch_file.write("N/A" + "\t")
                batch_file.write("3" + "\t")
                batch_file.write("N/A" + "\t")
                batch_file.write("N/A" + "\t")
                batch_file.write("N/A" + "\t")
                batch_file.write(cas_number + "\t")
                batch_file.write("Massbank EU" + "\n")

                for peak in peaks:
                    spectrum_string += peak + "\n"

                peaks = []
                spectrum_string += "END IONS\n"


                mgf_file.write(spectrum_string)

                if len(cas_number) < 2:
                    cas_number = "N/A"

                #batch_string = ""
                #batch_string += str(sys.argv[2]) + "\t"
                #batch_string += "Massbank:" + title + " " + compound_name + "\t"
                #batch_string += str(scan_number) + "\t"
                #batch_string += str(pepmass) + "\t"
                #batch_string += str(smiles) + "\t"
                #batch_string += str(inchi) + "\t"
                #batch_string += str(instrument) + "\t"
                #batch_string += str(exactmass) + "\t"
                #batch_string += str(cas_number) + "\t"
                #batch_string += str(adduct) + "\t"
                #batch_string += "\n"

                #batch_file.write(batch_string)

                scan_number += 1

                #Resetting variables
                peptide = "*..*"
                smiles = "N/A"
                inchi = "N/A"
                pepmass = ""
                title = ""
                instrument = ""
                compound_name = ""
                peaks = []
                retentiontime = ""
                ion_mode = ""
                exactmass = "0"
                cas_number = "N/A"
                adduct = "[M+H]"
                spectrum_level = 0


            #Resetting variables
            peptide = "*..*"
            smiles = "N/A"
            inchi = "N/A"
            pepmass = ""
            title = ""
            instrument = ""
            compound_name = ""
            peaks = []
            retentiontime = ""
            ion_mode = ""
            exactmass = "0"
            cas_number = "N/A"
            adduct = "[M+H]"
            spectrum_level = 0


        if line.find("ACCESSION") != -1:
            peptide = "*..*"
            peaks_start = 0
            title = line.replace("ACCESSION: ","").replace("//","").rstrip()
            #print line.replace("ACCESSION: ","")

        if line.find("CH$SMILES:") != -1:
            smiles = line[len("CH$SMILES: "):].rstrip()

        if line.find("CH$IUPAC: InChI=") != -1:
            inchi = line[len("CH$IUPAC: InChI="):].rstrip()

        if line.find("AC$MASS_SPECTROMETRY: ION_MODE") != -1:
            ion_mode = line[len("AC$MASS_SPECTROMETRY: ION_MODE "):].rstrip()

        if line.find("AC$INSTRUMENT_TYPE:") != -1:
            instrument = line[len("AC$INSTRUMENT_TYPE: "):].rstrip()

        if line.find("AC$CHROMATOGRAPHY: RETENTION_TIME ") != -1:
            retentiontime = line[len("AC$CHROMATOGRAPHY: RETENTION_TIME "):].rstrip()

        if line.find("MS$FOCUSED_ION: PRECURSOR_M/Z ") != -1:
            pepmass = line[len("MS$FOCUSED_ION: PRECURSOR_M/Z "):].rstrip()

        if line.find("MS$FOCUSED_ION: FULL_SCAN_FRAGMENT_ION_PEAK ") != -1:
            if len(pepmass) == 0:
                pepmass = line[len("MS$FOCUSED_ION: FULL_SCAN_FRAGMENT_ION_PEAK "):].rstrip()

        if line.find("CH$NAME: ") != -1:
            compound_name += line[len("CH$NAME: "):].rstrip() + "|"

        if line.find("CH$EXACT_MASS: ") != -1:
            exactmass = line[len("CH$EXACT_MASS: "):].rstrip()

        if line.find("CH$LINK: CAS ") != -1:
            cas_number = line[len("CH$LINK: CAS "):].rstrip()

        if line.find("MS$FOCUSED_ION: PRECURSOR_TYPE ") != -1:
            adduct = line[len("MS$FOCUSED_ION: PRECURSOR_TYPE "):].rstrip()

        if line.find("AC$MASS_SPECTROMETRY: MS_TYPE MS2") != -1:
            spectrum_level = 2

        if line.find("PK$PEAK") != -1:
            #print "PKPEAK"
            peaks_start = 1
            continue

        if (peaks_start == 1) and line.find("//") == -1:
            splits = line.split(" ")

            #print line
            peak_str = splits[2] +  " " + splits[3]
            if peak_str.find("int. rel.int") != -1:
                continue
            peaks.append(peak_str)
            #print peaks








if __name__ == "__main__":
    main()
