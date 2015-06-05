'''
Relatedness Estimator
Zheng Sun
CS124
'''

import numpy as np

def generate_genotypes_parents(maf, size):
    ''' generates genotypes for 1 parent
    args:
        maf = minor allele frequency // size = sequence size

    Genotype given by:
    0 = homogenous minor // 1 = homogenous major // 2 = heterogenous major
    '''
    p00 = maf * maf # probability of homogenous minor
    p01_10 = 2 * (1 - maf) * maf # probability of heterogenous
    p11 = (1 - maf) * (1 - maf) # probability of homogenous major

    genotypes = [] # store genotypes
    for i in range(0, size):
        p = np.random.uniform()
        if p < p00:
            genotypes.append(0)
        elif p < p00 + p01_10:
            genotypes.append(2)
        else:
            genotypes.append(1)
    return genotypes

def generate_genotypes_individuals(maf, size, n, related):
    ''' generates genotypes for n individuals
    args:
        maf = minor allele frequency // size = sequence size
        n = number of individuals' SNP to generate
        related = boolean of whether to generate related or unrelated individuals

    Genotype given by:
    0 = homogenous minor // 1 = homogenous major // 2 = heterogenous major
    '''

    # generate genotypes from each parent to pass down to child
    parent1 = generate_genotypes_parents(maf, size)
    parent2 = generate_genotypes_parents(maf, size)
    children = []
    for i in range(0, n):
        child = []
        for j in range(0, size):
            parent1_haplotype = parent1[j]
            parent2_haplotype = parent2[j]

            # if genotype is heterogenous, randomly select 0 or 1 as haplotype
            if parent1_haplotype == 2:
                parent1_haplotype = np.random.randint(2)
            if parent2_haplotype == 2:
                parent2_haplotype = np.random.randint(2)

            if parent1_haplotype == 0 and parent2_haplotype == 0: # homogenous minor
                child.append(0)
            elif parent1_haplotype == 0 and parent2_haplotype == 1: # heterogenous
                child.append(2)
            elif parent1_haplotype == 1 and parent2_haplotype == 0: # heterogenous
                child.append(2)
            elif parent1_haplotype == 1 and parent2_haplotype == 1: # homogenous major
                child.append(1)
        children.append(child)

        if not related: # generate new parents for each individual if not siblings
            parent1 = generate_genotypes_parents(maf, size)
            parent2 = generate_genotypes_parents(maf, size)
    return children

def baseline_comparison(a, b):
    counter = 0
    for i in range(len(a)):
        if a[i] == b[i]:
            counter += 1
        else:
            counter -= 1
    return (counter > 0)

def run_baseline(maf, size, n):
    correct = 0
    for i in range(0, n):
        if np.random.randint(10) >= 7:
            compared_SNPs = generate_genotypes_individuals(maf, size, 2, True)
            if baseline_comparison(compared_SNPs[0], compared_SNPs[1]) == True:
                correct += 1
        else:
            compared_SNPs = generate_genotypes_individuals(maf, size, 2, False)
            if baseline_comparison(compared_SNPs[0], compared_SNPs[1]) == False:
                correct += 1
    return correct / float(n)

def improved_comparison(a, b, p):
    unrelated = [[0 for x in range(3)] for x in range(3)]
    related = [[0 for x in range(3)] for x in range(3)]

    '''unrelated[0][0] = (1-p) * (1-p) * (1-p) * (1-p)
    unrelated[0][1] = (1-p) * (1-p) * p * p
    unrelated[0][2] = 2 * (1-p) * (1-p) * (1-p) * p
    unrelated[1][0] = (1-p) * (1-p) * p * p
    unrelated[1][1] = p * p * p * p
    unrelated[1][2] = 2 * (1-p) * p * p * p
    unrelated[2][0] = 2 * (1-p) * (1-p) * (1-p) * p
    unrelated[2][1] = 2 * (1-p) * p * p * p
    unrelated[2][2] = 4 * (1-p) * (1-p) * p * p

    related[0][0] = unrelated[0][0] * 1 + unrelated[0][2] * 0.25 + unrelated[2][0] * 0.25 + unrelated[2][2] * 0.0625
    related[0][1] = unrelated[2][2] * 0.0625
    related[0][2] = unrelated[0][2] * 0.25 + unrelated[2][0] * 0.25 + unrelated[2][2] + 0.125
    related[1][0] = unrelated[2][2] * 0.0625
    related[1][1] = unrelated[1][1] * 1 + unrelated[1][2] * 0.25 + unrelated[2][1] * 0.25 + unrelated[2][2] * 0.0625
    related[1][2] = unrelated[1][2] * 0.25 + unrelated[2][1] * 0.25 + unrelated[2][2] * 0.125
    related[2][0] = unrelated[0][2] * 0.25 + unrelated[2][0] * 0.25 + unrelated[2][2] * 0.125
    related[2][1] = unrelated[1][2] * 0.25 + unrelated[2][1] * 0.25 + unrelated[2][2] * 0.125
    related[2][2] = unrelated[0][1] * 1 + unrelated[1][0] * 1 + unrelated[0][2] * 0.25 + unrelated[2][0] * 0.25 + unrelated[1][2] * 0.25 + unrelated[2][1] * 0.25 + unrelated[2][2] * 0.25
    '''

    unrelated[0][0] = (1-p) * (1-p) * (1-p) * (1-p)
    unrelated[0][1] = p * p * (1-p) * (1-p)
    unrelated[0][2] = 2 * p * (1-p) * (1-p) * (1-p)
    unrelated[1][0] = p * p * (1-p) * (1-p)
    unrelated[1][1] = p * p * p * p
    unrelated[1][2] = 2 * p * p * p * (1-p)
    unrelated[2][0] = 2 * p * (1-p) * (1-p) * (1-p)
    unrelated[2][1] = 2 * p * p * p * (1-p)
    unrelated[2][2] = 4 * p * p * (1-p) * (1-p)

    related[0][0] = p * p * (1-p) * (1-p) * (.25) + p * (1-p) * (1-p) * (1-p) + (1-p) * (1-p) * (1-p) * (1-p)
    related[0][1] = p * p * (1-p) * (1-p) * (.25)
    related[0][2] = p * p * (1-p) * (1-p) * (.5) + p * (1-p) * (1-p) * (1-p)
    related[1][0] = p * p * (1-p) * (1-p) * (0.25)
    related[1][1] = p * p * p * p + p * p * p * (1-p) + p * p * (1-p) * (1-p) * (.25)
    related[1][2] = p * p * p * (1-p) + p * p * (1-p) * (1-p) * (0.5)
    related[2][0] = p * p * (1-p) * (1-p) * (0.5) + p * (1-p) * (1-p) * (1-p)
    related[2][1] = p * p * p * (1-p) + p * p * (1-p) * (1-p) * (0.5)
    related[2][2] = p * p * p * (1-p) + p * p * (1-p) * (1-p) * 3 + p * (1-p) * (1-p) * (1-p)

    num_unrelated = 0
    num_related = 0

    for i in range(len(a)):
        if unrelated[a[i]][b[i]] > related[a[i]][b[i]]:
            num_unrelated += 1
        #elif unrelated[a[i]][b[i]] <= related[a[i]][b[i]]:
        else:
            num_related += 1
        
    return num_related >= num_unrelated

def run_improved(maf, size, n):
    correct = 0
    for i in range(0, n):
        if np.random.randint(10) >= 7:
            compared_SNPs = generate_genotypes_individuals(maf, size, 2, True)
            if improved_comparison(compared_SNPs[0], compared_SNPs[1], maf) == True:
                correct += 1
        else:
            compared_SNPs = generate_genotypes_individuals(maf, size, 2, False)
            if improved_comparison(compared_SNPs[0], compared_SNPs[1], maf) == False:
                correct += 1
    return correct / float(n)

print run_baseline(0.4, 100, 1000)
print run_improved(0.4, 100, 1000)