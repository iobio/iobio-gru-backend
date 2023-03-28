#!/usr/bin/env python3
# Given a human gene symbol and species ids, returns ensemble gene ids
# for orthologs corresponding to provided species ids, if they exist.
# If no species ids provided, just pulls back zebrafish, mouse, and human.
# If not orthologs could be found, returns None
# SJG for Gabor Marth lab Mar2023

import requests
import time


# Returns a list of human gene IDs from NCBI corresponding to a given gene name
def get_ncbi_gene_ids(gene_name):
    eutils_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    iobio_eutils_key = "2ce5a212af98a07c6e770d1e95b99a2fef09"

    params = {'db': 'gene',
              'retmode': 'json',
              'term': '(9606[Taxonomy ID] AND ' + gene_name + '[Pref])',
              'api_key': iobio_eutils_key
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
def get_ortho_ids(gene_name, gene_ids, species_ids):
    odb_api_url = "https://data.orthodb.org/current/"
    zebrafish_eutils_id = "7955"
    mouse_eutils_id = "10090"
    human_eutils_id = "9606"
    vertebrata_tax_id = 7742  # Highest level of clade that human, mouse, zebrafish converge at in OrthoDB

    if species_ids is None:
        species_ids = [zebrafish_eutils_id, mouse_eutils_id, human_eutils_id]

    cluster_ids = []
    for gene_id in gene_ids:
        params = {'query': gene_name + " " + gene_id,
                  'level': vertebrata_tax_id
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
        if not payload:
            print("Incorrect response from orthoDB. Exiting...")
            return None

        data = payload.get("data")
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
    if not payload:
        print("Incorrect response from uniprot. Exiting...")
        return None
    results = payload.get("results")
    if not results:
        print("Empty results from uniprot. Exiting...")
        return []
    ensembl_ids = [result["to"][:result["to"].index('.')] for result in payload["results"]]
    return ensembl_ids


def print_orthologs(gene_name, species_ids=None):
    gene_ids = get_ncbi_gene_ids(gene_name)
    if not gene_ids:
        return
    uniprot_ids = get_ortho_ids(gene_name, gene_ids, species_ids)
    if not uniprot_ids:
        return
    ensembl_ids = get_ensembl_ids(uniprot_ids)
    if not ensembl_ids:
        return
    print(ensembl_ids)


if __name__ == '__main__':
    sample_gene = 'APC'
    print_orthologs(sample_gene)
