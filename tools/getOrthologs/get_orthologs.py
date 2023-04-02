#!/usr/bin/env python3
# Given a human gene symbol and species ids, returns ensembl gene ids
# for orthologs corresponding to provided species ids, if they exist.
# If no orthologs could be found, returns None.
# Required arguments: gene name, species ID list (NCBI assigned), taxon level ID (NCBI assigned), and eutils key.
# NOTE: the taxon level ID corresponds to the narrowest evolutionary convergence level that
# encompasses the provided species. When using Chordata for mouse/zebrafish/human, however, orthoDB failed.
# Thus using Verterbrata for that specific list. Any other level should be tested directly with orthoDB
# before being provided as an argument here.
# SJG for Gabor Marth lab Mar2023

import requests
import time
import sys
import getopt


def get_args(argv):
    usage = 'get_orthologs.py -g <geneName> -k <eutilsKey> -s <speciesIds> -t <taxonLevelId>'
    try:
        opts, args = getopt.getopt(argv, "hg:k:s:t:", ["geneName=", "eutilsKey=", "speciesIds=", "taxonLevelId="])
    except getopt.GetoptError:
        print(usage, file=sys.stderr)
        sys.exit(2)

    the_args = {}
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-g", "--geneName"):
            the_args['gene_name'] = arg
        elif opt in ("-k", "--eutilsKey"):
            the_args['eutils_key'] = arg
        elif opt in ("-s", "--speciesIds"):
            the_args['species_ids'] = arg
        elif opt in ("-t", "--taxonLevelId"):
            the_args['taxon_level_id'] = arg

    return the_args


# Returns a list of human gene IDs from NCBI corresponding to a given gene name
def get_ncbi_gene_ids(the_args):
    gene_name = the_args['gene_name']
    eutils_key = the_args['eutils_key']

    eutils_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    params = {'db': 'gene',
              'retmode': 'json',
              'term': '(9606[Taxonomy ID] AND ' + gene_name + '[Pref])',
              'api_key': eutils_key
              }
    r = requests.get(eutils_url, params=params)

    payload = r.json()
    result = payload.get("esearchresult")
    if payload is None or not result:
        print("Incorrect response from eutils. Exiting...")
        return None

    id_list = result.get("idlist")
    if not id_list:
        print("Empty response from eutils. Exiting...")

    return id_list


# Returns UNIPROT IDs from ortho-db corresponding to given NCBI gene IDs
# If no response from UNIPROT, returns None
# https://www.ezlab.org/orthodb_userguide.html#api for more info
def get_ortho_ids(the_args, gene_ids):
    gene_name = the_args['gene_name']
    species_ids = the_args['species_ids']
    taxon_level_id = the_args['taxon_level_id']

    odb_api_url = "https://data.orthodb.org/current/"

    cluster_ids = []
    for gene_id in gene_ids:
        params = {'query': gene_name + " " + gene_id,
                  'level': taxon_level_id
                  }
        r = requests.get(odb_api_url + "search", params=params)
        payload = r.json()
        data = payload.get("data")
        if data:
            cluster_ids.extend(data)

    if not cluster_ids:
        print("Could not find cluster for gene " + gene_name + ". Exiting...")
        return None

    uniprot_ids = []
    for cluster_id in cluster_ids:
        params = {'id': cluster_id}
        r = requests.get(odb_api_url + "orthologs", params=params)
        payload = r.json()
        data = payload["data"]
        if data:
            gene_lists = [ortholog["genes"] for ortholog in data if
                          (ortholog["organism"]["id"][:ortholog["organism"]["id"].index('_')]) in species_ids]
            for gene_list in gene_lists:
                for gene in gene_list:
                    if gene.get("uniprot"):
                        uniprot_ids.append(gene["uniprot"]["id"])
        if data is None:
            print("Incorrect response from orthoDB. Exiting...")
            return None

    uniq_ids = list(set(uniprot_ids))
    return uniq_ids


# Returns list of ensembl IDs corresponding to the given UNIPROT IDs
# If no response from Uniprot, returns None
# https://www.uniprot.org/help/id_mapping for more info
def get_ensembl_ids(uniprot_ids):
    uniprot_api_url = "https://rest.uniprot.org/idmapping/"
    params = {'ids': ','.join(uniprot_ids),
              'from': 'UniProtKB_AC-ID',
              'to': 'Ensembl'}
    r = requests.post(uniprot_api_url + "run", params=params)
    payload = r.json()
    if not payload:
        print("Incorrect response from uniprot. Exiting...")
        return None

    job_id = payload.get("jobId")
    if not job_id:
        print("Incorrect response from uniprot. Exiting...")
        return None

    # Poll until results ready (or we've tried 30 times)
    uniprot_working = True
    count = 0
    r = requests.get(uniprot_api_url + "status/" + str(job_id))
    while uniprot_working and count < 30:
        payload = r.json()
        if payload.get("results"):
            uniprot_working = False
        elif payload.get("jobStatus"):
            uniprot_working = payload["jobStatus"] == "FINISHED"
        time.sleep(1)
        count += 1

    r = requests.get(uniprot_api_url + "results/" + job_id)
    payload = r.json()
    ensembl_ids = [result["to"][:result["to"].index('.')] for result in payload["results"]]
    return ensembl_ids


def print_orthologs(the_args):
    gene_ids = get_ncbi_gene_ids(the_args)
    if not gene_ids:
        sys.exit(1)
    uniprot_ids = get_ortho_ids(the_args, gene_ids)
    if not uniprot_ids:
        sys.exit(1)
    ensembl_ids = get_ensembl_ids(uniprot_ids)
    if not ensembl_ids:
        sys.exit(1)
    print(ensembl_ids)


if __name__ == '__main__':
    args = get_args(sys.argv[1:])
    print_orthologs(args)
