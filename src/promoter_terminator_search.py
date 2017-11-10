import urllib.request as urllib
import shutil
import zipfile
import os

class PromoterTerminatorSearch:
    def __init__(self):
        #hardcoded url
        self.request_url = open('../bin/restriction_promoter_url.txt', 'r').read()

    def get_data(self, url, output_filename):
        with urllib.urlopen(url) as response, open(output_filename, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)

if __name__ == '__main__':
    p = PromoterTerminatorSearch()
    url = 'http://pallab.serc.iisc.ernet.in/gester/05712904/RawData.zip'
    p.get_data(url, 'RawData.zip')
    os.system('unzip RawData.zip')
