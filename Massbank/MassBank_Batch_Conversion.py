#!/usr/bin/python


import sys
import getopt
import os



def usage():
    print "<input folder>"




def main():
    usage()

    input_folder_path = sys.argv[1]
    
    files = os.listdir(input_folder_path)
    
    reported_fields = []
    reported_fields.append("NAME")
    reported_fields.append("SMILES")
    reported_fields.append("INCHI")
    reported_fields.append("SOURCEINSTRUMENT")
    reported_fields.append("ORGANISM")
    reported_fields.append("PEPMASS")
    
    for input_file in files:
        fileName, fileExtension = os.path.splitext(input_file)
        if fileExtension.find("mgf") == -1:
            continue
        
        output_temp = ""

        spectra_count = 0
        charge = 0
        precursor_mz = 0.0
        for line in open(input_file, "r"):
            if line.find("SCANS=") != -1:
                revised_file.write("SCANS=" + str(spectra_count + 1) + "\n")
                continue
            revised_file.write(line.rstrip() + "\n")
            if line.find("CHARGE=") != -1:
                charge = int(line.rstrip()[7:])
            if line.find("PEPMASS") != -1:
                precursor_mz = float(line.rstrip()[8:])
            if line.find("NAME=") != -1:
                output_temp += line.rstrip()[5:] + "\t" 
            if line.find("NAME=") != -1:
                output_temp += line.rstrip()[5:] + "\t" 
            if line.find("NAME=") != -1:
                output_temp += line.rstrip()[5:] + "\t" 
                
            if line.find("END IONS") != -1:
                spectra_count += 1
                print input_file + "\t" + str(charge) + "\t" + str(precursor_mz) + "\t" + str(spectra_count) + "\t" + output_temp
            
        
    #print input_file + "\t" + str(spectra_count)
    
    

if __name__ == "__main__":
    main()
