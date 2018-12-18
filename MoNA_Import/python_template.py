#!/usr/bin/python


import sys
import getopt
import os
import requests
import json
import ming_parallel_library

def usage():
    print "<output json folder>"

def get_total_spectra():
    url = "http://mona.fiehnlab.ucdavis.edu/rest/spectra/searchCount"
    response = requests.get(url)
    json_obj = json.loads(response.text)
    return int(json_obj["count"])

def main():
    output_json_folder = sys.argv[1]

    entries = 1000

    total_spectra = get_total_spectra()

    input_data = []

    for i in range(total_spectra/entries + 1):
        offset = i * entries
        input_data.append({"entries": entries, "offset": offset, "json_folder": output_json_folder})

        #print "Grabbing offset: " + str(offset)
        #spectra_data = pulldown_mona(entries, offset)
        #for spectrum in spectra_data:
        #    print spectrum["id"]

    results = ming_parallel_library.run_parallel_job(pulldown_mona_parallel_wrapper, input_data, 10)

    #while(True):
    #    print "Grabbing offset: " + str(offset)
    #    spectra_data = pulldown_mona(entries, offset)
    #    #print spectra_data
    #    #print len(spectra_data)
    #    for spectrum in spectra_data:
    #        print spectrum["id"]
    #    offset += len(spectra_data)
    #    if len(spectra_data) == 0:
    #        break


def pulldown_mona_parallel_wrapper(parameters):
    entries = parameters["entries"]
    offset = parameters["offset"]

    json_obj = pulldown_mona(entries, offset)

    output_folder = parameters["json_folder"]
    output_filename = os.path.join(output_folder, str(entries) + "_" + str(offset) + ".json")
    output_file = open(output_filename, "w")
    output_file.write(json.dumps(json_obj))

    return 0

def pulldown_mona(entries, offset):
    print "Getting Data Offset: " + str(offset)
    payload = {"max": entries, "offset": offset}
    url = "http://mona.fiehnlab.ucdavis.edu/rest/spectra/search?max=" + str(entries) + "&offset=" + str(offset)
    response = requests.get(url, params='{"compound":{},"metadata":[],"tags":[],"offset":' + str(offset) + '}')
    json_obj = json.loads(response.text)
    print "Got: " + str(len(json_obj)) + " spectra"

    return json_obj



if __name__ == "__main__":
    main()
