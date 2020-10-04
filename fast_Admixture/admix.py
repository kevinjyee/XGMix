"""
Admix++
- Get a faster pipe for getting admixture simulation.
- Follow an approximate method to get this done.
- Allow weights for different ancestries - 
    seen this use case in our xgmix, lainet work.
- Get numpy files directly and be mindful of storage requirements
- sparse for vcf, rle for ancestries could be future steps

"""

import allel
import sys

from utils import build_founders, create_dataset, write_output

def simulate(reference, sample_map, genetic_map, out_root,
             num_samples_per_gen, gens_to_ret,
             sample_weights=None,random_seed=42):

    """

    out_root is modified to infer chromosome number, gen number. 

    """

    print("Reading vcf data")
    vcf_data = allel.read_vcf(reference)
    print("Building founders")
    founders, founders_weight= build_founders(vcf_data, 
                                               sample_map, 
                                               genetic_map, 
                                               sample_weights)
    chm = list(set(vcf_data["variants/CHROM"]))[0]
    print("Simulating...")
    dataset = create_dataset(founders,founders_weight,num_samples_per_gen,gens_to_ret,random_seed)
    print("Writing output")
    out_root = out_root + "/chm"+str(chm)
    write_output(out_root,dataset)

if __name__ == "__main__":

    """
    Sample command:
    python admix.py /home/arvindsk/datagen/world_wide_references/ref_final_beagle_phased_1kg_hgdp_sgdp_chr22.vcf /home/wknd37/Admixture/generated_data/6_even_anc_t2/chm22/sample_maps/train.map /home/database/maps/rfmix/allchrs.b37.gmap /home/arvindsk/datagen/fast_admix/chm22/train/ 8 400

    """

    reference = sys.argv[1]
    sample_map = sys.argv[2]
    genetic_map = sys.argv[3]
    out_root = sys.argv[4]

    max_gen = int(sys.argv[5])
    # returns 2 to max_gen.
    # by default we write gen 0.
    gens_to_ret = range(2,max_gen+1)

    num_samples_per_gen = int(sys.argv[6])

    random_seed=42
    sample_weights=None

    if len(sys.argv) >= 8:
        random_seed = int(sys.argv[7])

    if len(sys.argv) >= 9:
        sample_weights = sys.argv[8]

    simulate(reference, sample_map, genetic_map, out_root,
             num_samples_per_gen, gens_to_ret,
             sample_weights,random_seed)
