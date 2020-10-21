"""
utils.py contains helper function to deal with creating founders,
processing genetic map and processing sample maps.
"""

import numpy as np
import allel
import pandas as pd
import os

from .person import Person, create_new

def get_chm_lengths(genetic_map):

    """
    get chromosome length in morgans from genetic map.
    """

    chromosome_lengths_morgans = {}
    with open(genetic_map, "r") as f:
        for line in f.readlines():
            a,b,c = line.split("\t")
            chromosome_lengths_morgans[a] = float(c)/100 # cM to M
            
    # print(chromosome_lengths_morgans)
    return chromosome_lengths_morgans

def get_sample_map_data(sample_map, vcf_data, sample_weights = None):
    
    """
    Inputs:
    sample_map: tab delimited file with sample, population and no header.
    vcf_data: allel.read_vcf output. It is the reference vcf file information.
    
    Returns:
    sample_map_data: dataframe with sample, population, population_code and index in vcf_data referecnce.
    
    """
    
    # reading sample map
    sample_map_data = pd.read_csv(sample_map,delimiter="\t",header=None)
    sample_map_data.columns = ["sample","population"]

    # creating ancestry map into integers from strings
    ancestry_map = {}
    curr = 0
    for i in sample_map_data["population"]:
        if i in ancestry_map.keys():
            continue
        else:
            ancestry_map[i] = curr
            curr += 1
    print("Ancestry map",ancestry_map)
    sample_map_data["population_code"] = np.vectorize(ancestry_map.get)(sample_map_data["population"])

    # getting index of samples in the reference files

    b = vcf_data["samples"]
    a = np.array(list(sample_map_data["sample"]))

    sorter = np.argsort(b)
    indices = sorter[np.searchsorted(b, a, sorter=sorter)]
    sample_map_data["index_in_reference"] = indices
    
    if sample_weights is not None:
        sample_weights_df = pd.read_csv(sample_weights,delimiter="\t",header=None)
        sample_weights_df.columns = ["sample","sample_weight"]
        pd.merge(sample_map_data, sample_weights_df, on='sample')
        
    else:
        sample_map_data["sample_weight"] = [1.0/len(sample_map_data)]*len(sample_map_data)
    
    return sample_map_data

def build_founders(vcf_data, sample_map, genetic_map, sample_weights = None):
    
    """
    Returns founders - a list of Person datatype.
    founders_weight - a list with a weight for each sample in founders
    
    """
    
    
    sample_map_data = get_sample_map_data(sample_map, vcf_data)
    
    chromosome_lengths_morgans = get_chm_lengths(genetic_map)

    chm = list(set(vcf_data["variants/CHROM"]))
    if len(chm) == 1:
        chm = chm[0]
    else:
        raise Exception("Multiple chromosomes in this file!!!")


    chm_length_morgans = chromosome_lengths_morgans[chm]
    chm_length_snps = vcf_data["calldata/GT"].shape[0]


    # building founders
    founders = []

    for i in sample_map_data.iterrows():

        # first get the index of this sample in the vcf_data.
        # if not there, skip and print to log.

        index = i[1]["index_in_reference"]

        name = i[1]["sample"]

        # when creating maternal, paternal make sure it has same keys


        maternal = {}
        paternal = {}

        # let us use the first for maternal in the vcf file...
        maternal["snps"] = vcf_data["calldata/GT"][:,index,0]
        paternal["snps"] = vcf_data["calldata/GT"][:,index,1]

        # single ancestry assumption.
        maternal["anc"] = np.array([i[1]["population_code"]]*chm_length_snps)
        paternal["anc"] = np.array([i[1]["population_code"]]*chm_length_snps)


        # any more info like coordinates, prs can be added here.

        p = Person(chm,chm_length_morgans,chm_length_snps,maternal,paternal,name)

        founders.append(p)
        
    founders_weight = sample_map_data["sample_weight"]
    
    return founders,founders_weight

def create_dataset(founders,founders_weight,num_samples_per_gen,gens_to_ret,random_seed=42):
    
    np.random.seed(random_seed)
    
    max_gen = max(gens_to_ret)
    
    overall = {}
    overall[0] = founders

    for gen in range(1,max_gen+1):
        print("Simulating generation ",gen)
        this_gen_people = []


        for i in range(num_samples_per_gen):
            # select any 2 parents from prev. gen.
            # if 1st generation, choose from founders
            if gen == 1:
                id1,id2 = np.random.choice(len(founders),size=2,replace=False, p=founders_weight)
            # else choose from num_samples_per_gen
            else:
                id1,id2 = np.random.choice(num_samples_per_gen,size=2,replace=False)
            p1 = overall[gen-1][id1]
            p2 = overall[gen-1][id2]

            adm = create_new(p1,p2)
            this_gen_people.append(adm)

        overall[gen] = this_gen_people

        if gen-1 not in gens_to_ret and gen-1 !=0:
            del overall[gen-1]

    return overall

def write_output(root,dataset):
    
    """
    creates output numpy files in root directory - under that, we will have gen_{}
    folders and then we have the npy files.
    Make sure root has chm number in it.
    
    """

    if not os.path.isdir(root):
        os.makedirs(root)

    for gen in dataset.keys():
        
        print("Writing generation: {}".format(gen))

        if gen not in dataset.keys():
            print("Did not simulate gen {}".format(gen))
            continue

        # make directory
        gen_dir = os.path.join(root, "gen_{}".format(gen))
        if not os.path.isdir(gen_dir):
            os.makedirs(gen_dir)

        # create npy files.

        mat_snps = np.vstack([i.maternal["snps"] for i in dataset[gen]])
        pat_snps = np.vstack([i.paternal["snps"] for i in dataset[gen]])
        snps = np.vstack([mat_snps,pat_snps])
        np.save(gen_dir+"/mat_vcf_2d.npy",snps)

        # create map files.

        mat_anc = np.vstack([i.maternal["anc"] for i in dataset[gen]])
        pat_anc = np.vstack([i.paternal["anc"] for i in dataset[gen]])
        anc = np.vstack([mat_anc,pat_anc])
        np.save(gen_dir+"/mat_map.npy",anc)
     
    