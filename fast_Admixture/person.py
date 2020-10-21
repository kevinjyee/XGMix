"""
- Person class definition.

A person with maternal and paternal sequence information.
Has an admix method to create next generation admixed individual.

- create_new function

Takes 2 Person datatypes and returns an admixed Person

"""

import numpy as np

class Person():
    
    def __init__(self,chm,chm_length_morgans,chm_length_snps,maternal,paternal,name=None):
          
        """
        Inputs:
        chm: chm number.
        chm_length_morgans: chromosome length in morgans
                found using the genetic map.
        chm_length_snps: chromosome length in avaliable snps
                based on reference file.
        
        maternal, paternal are dictionaries with the keys: snps, anc, etc...
            this gives us easy extension to co-ordinates, prs, etc...
        
        """

        # christen the admixed individual. :P
        self.name = name if name is not None else \
                "admixed"+str(int(np.random.rand()*1e6))
        
        # chm related information
        self.chm = chm
        self.chm_length_morgans = chm_length_morgans
        self.chm_length_snps = chm_length_snps
        
        assert maternal.keys() == paternal.keys(), "Key mismatch!!!"
        
        self.order = sorted(maternal.keys())
        
        self.maternal = maternal
        self.paternal = paternal    

        # assert all sequences have same length
        
    def admix(self):

        """
        create an admixed haploid from the paternal and maternal sequences.
        
        Returns:
        haploid_returns: dict with same keys as self.maternal and self.paternal

        """

        num_crossovers = int(np.random.poisson(self.chm_length_morgans))
        
        # debug...
        # num_crossovers=3
        
        #print(num_crossovers)
        
        haploid_returns = {}
        for key in self.order:
            haploid_returns[key] = np.zeros_like(self.maternal[key])
        
        # edge case of no breaking points.
        if num_crossovers == 0:
                
            haploid_returns = {}
            select = self.maternal if np.random.rand()>=0.5 else self.paternal
            for key in self.order:
                
                haploid_returns[key] = select[key]
                
        
        # now, we have to choose crossover points and return 
        # the appropriate admixed sequence.
        
        # the crossover points must be chosen based on the 
        # genetic map and not the snps we have.
        # TODO: gotta make sure to pick the closest point in the snps we have.
        # For now, just uniform over snp length.
        
        else:
            
            breakpoints = np.random.choice(np.arange(1,self.chm_length_snps), 
                                           size=num_crossovers, 
                                           replace=False, 
                                           p=None)
            breakpoints = np.sort(breakpoints)
            
            breakpoints = np.concatenate(([0],breakpoints,[self.chm_length_snps]))
            
            #print(breakpoints)
            # select paternal or maternal randomly and apply crossovers.
            choice = np.random.rand()>=0.5
            select = self.maternal if choice else self.paternal
            cross = self.paternal if choice else self.maternal
            for key in self.order:
                haploid_returns[key] = select[key].copy()
                for i in range(1,len(breakpoints)-1,2):
                    begin = breakpoints[i]
                    end = breakpoints[i+1]
                    haploid_returns[key][begin:end] = cross[key][begin:end]
    
        return haploid_returns


def create_new(p1,p2):
    
    maternal = p1.admix()
    paternal = p2.admix()
    
    assert p1.chm == p2.chm, "Wrong chromosomes being passed!!!"
    
    chm = p1.chm
    chm_length_morgans = p1.chm_length_morgans
    chm_length_snps = p1.chm_length_snps

    return Person(chm,chm_length_morgans,chm_length_snps,maternal,paternal)