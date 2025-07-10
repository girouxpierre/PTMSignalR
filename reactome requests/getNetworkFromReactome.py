import json
import os
from utils.parsersExport import parseExportComplex, parseExportTransport, parseExportPhospho, parseExportExpression, parseExportCatalysisPrecedes
from utils.parsersExportStateOf import parseExportStateOf
from utils.requests import requestComplexes, requestTransport, requestStateOf, requestPhospho, requestExpression, requestCatalysisPrecedes
from utils.functions import convertUniprotToGeneName, printLogo, mergingAllSif

def main():
    '''

    - set up config.json : set complexeToIgnore.txt Path and metabolistToRM.txt path
    OR
    - set your own threshold for complex size by launching getComplexeSize.txt
    - set your own treshold for metabolic links by directly modyfying request in configRequest.Use countMetabolitesRelationship.py to get your own list of ID
    '''

    print(f' ==========     reading Config   ==========')

    with open('src/config.json', 'r') as config_file:
        config = json.load(config_file)

    uri = config['uri']
    user = config['user']
    password = config['password']
    complexesToIgnorePath = config["complexesToIgnorePath"]
    print("Launching with ", complexesToIgnorePath, "as complexes To remove\n")
    HUGOBiomartNamePath = config["HugoRefFile"]
    jsonConfigRequestPath = config["jsonConfigRequestPath"]
    metaboToIgnored = config["metaboToIgnore50"]
    print("Launching with ", jsonConfigRequestPath, "as requestFile with metabolites to ignore\n")

    with open(jsonConfigRequestPath, 'r') as config_file2:
        configRequests = json.load(config_file2)


    cypPhospho_genePin = configRequests['cypPhospho_genePin']
    cypPhospho_genePout = configRequests['cypPhospho_genePout']

    # control Phospho of 

    if not(os.path.isfile('out/sif/controlPhospho.sif')):
        print(f' ==========     getting ControlPhosphoOf.sif   ==========')
        requestPhospho(uri, user, password, cypPhospho_cataSolo, cypPhospho_cataComp, cypPhospho_genePin, cypPhospho_genePout, 'tmp/phospho')
        parseExportPhospho('tmp/phospho_CataSolo.tsv', 'tmp/phospho_CataComp.tsv', 'tmp/stateOf_GeneReg.tsv', 'tmp/phospho_GenePin.tsv', 'tmp/phospho_GenePout.tsv', 'tmp/controlPhopshoUniprot.sif')
        convertUniprotToGeneName(HUGOBiomartNamePath, 'tmp/controlPhopshoUniprot.sif', 'out/sif/controlPhospho.sif')

    # merge all sif
    
    mergingAllSif('out/sif/', 'extriTF/TF_target.sif', 'data/manuallyAddedInteractions.sif','out/final/network.sif') # merge all sif and integrates TF-Target (+ manual interaction ?) 

# Execute
if __name__ == "__main__":
    printLogo('assets/logo.txt')
    main()
