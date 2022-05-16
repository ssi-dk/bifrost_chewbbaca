#from urllib import HTTPError
import urllib
import json
import os
from urllib.error import HTTPError
import urllib.request
import argparse
import subprocess
import html
import bs4
import re
from multiprocessing import Pool
parser = argparse.ArgumentParser()
parser.add_argument("--address", "-a", default = "https://enterobase.warwick.ac.uk/schemes/Yersinia.cgMLSTv1/", help = "URL containing enterobase fasta links")
parser.add_argument("--output_folder", "-o", default = "/srv/data/MPV/LEBC/test_download", help = "scheme output folder")
args = parser.parse_args()
output_folder = args.output_folder
#YERSINIA_ADDRESS = 'https://enterobase.warwick.ac.uk/schemes/Yersinia.cgMLSTv1/'
#DATABASE = 'yersinia'
#scheme = 'cgMLST_v2'
API_TOKEN = "eyJhbGciOiJIUzI1NiIsImV4cCI6MTY2MjUyMTI5MiwiaWF0IjoxNjQ2NzUzMjkyfQ.eyJ1c2VybmFtZSI6Ik1hb3MiLCJjaXR5IjoiQ29wZW5oYWdlbiIsImNvbmZpcm1lZCI6MSwiYWxsb3dlZF9zY2hlbWVzIjoiY2dNTFNUX3YyIiwiZmlyc3RuYW1lIjoiTWFyayIsImFwaV9hY2Nlc3Nfc2VudGVyaWNhIjoiVHJ1ZSIsImNvdW50cnkiOiJEZW5tYXJrIiwiaWQiOjEzMzAsImFkbWluaXN0cmF0b3IiOm51bGwsImVtYWlsIjoibWFvc0Bzc2kuZGsiLCJkZXBhcnRtZW50IjoiRm9vZGJvcm5lIHBhdGhvZ2VuIHN1cnZlaWxsYW5jZSIsInZpZXdfc3BlY2llcyI6IlRydWUiLCJsYXN0bmFtZSI6Ik9lc3Rlcmx1bmQiLCJhY3RpdmUiOm51bGwsInVwbG9hZF9yZWFkcyI6IlRydWUiLCJpbnN0aXR1dGlvbiI6IlN0YXRlbnMgU2VydW0gSW5zdGl0dXQifQ.bIrgc4P90xBAFlBCoPVOHQ-Ui0PqYWKZmGHZJdmLGMk"

if not os.path.isdir(output_folder):
    os.mkdir(output_folder)

def get_html(address):
    curl_cmd = f"curl -k -X GET --header 'Accept: application/json' --user '{API_TOKEN}:' {address}"
    html_output = subprocess.check_output(f"curl -k -X GET --header 'Accept: application/json' --user '{API_TOKEN}:' {address}", shell=True)
    soup = bs4.BeautifulSoup(html_output, features='html.parser')
    links_html = soup.find_all("a")
    links = [os.path.join(address, i['href']) for i in links_html if re.match(".*fasta.gz", i['href'])]
    return links
    
def download_fastas(url):
    return_val = None
    try:
        u = urllib.request.urlopen(url)
    except urllib.error.URLError as error:
        return_val = error
    output_path = os.path.join(output_folder, url.split("/")[-1])
    with open(output_path, "wb+") as output:
        output.write(u.read())
    return return_val==None

def repeated_download(url):
    N_retries = 2
    while N_retries > 0:
        download_bool = download_fastas(url)
        if not download_bool:
            print(f"dowloading {url} failed, {N_entries} attempts left")
        else:
            break
        
#full_address = f"{API_ADDRESS}{DATABASE}/{scheme}/loci?limit=0&offset=0"
#cmd = f"curl -k -X GET --header 'Accept: application/json' --user '{API_TOKEN}:' {full_address}"
#process  = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, env=os.environ, encoding='utf-8')
#process_out, process_err = process.communicate()
#print(get_N_loci(), "output")
#print(process_out, "process_out")
print('getting fasta links')
fasta_links = get_html(args.address)
pool = Pool(processes=4)
print('downloading fastas with 4 cores')
pool.map(download_fastas, fasta_links)
#for i in fasta_links:
#    download_fastas(i)
#download_fastas(fasta_links[0])
#print(fasta_links[0], len(fasta_links))
