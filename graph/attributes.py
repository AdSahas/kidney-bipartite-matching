import random
random.seed(42)  # for reproducibility and statistical consistency across runs.

# distrubutions: modelled after US/European blood types
BLOOD_DISTRIBUTION = {'O': 0.44, 'A': 0.42, 'B': 0.10, 'AB': 0.04}
RH_DISTRIBUTION = {'+': 0.85, '-': 0.15}

def generateBloodType(n):
    bloodType = random.choices(
        list(BLOOD_DISTRIBUTION.keys()), 
        weights=list(BLOOD_DISTRIBUTION.values()),
        k = n
        )
    return bloodType

def generateRhFactor(n):
    rhFactor = random.choices(
        list(RH_DISTRIBUTION.keys()), 
        weights=list(RH_DISTRIBUTION.values()),
        k = n
    )
    return rhFactor


# HLA and PRA for transplant compatibility
# https://www.ebi.ac.uk/ipd/imgt/hla/

# Robinson J, et al. "The IPD and IMGT/HLA database: allele variant databases." 
# Nucleic Acids Research, 2015. https://www.ebi.ac.uk/ipd/imgt/hla/
HLA_A_ALLELES = [
    'A1', 'A2', 'A3', 'A11', 'A23', 'A24', 'A25', 'A26',
    'A29', 'A30', 'A31', 'A32', 'A33', 'A34', 'A36', 'A43',
    'A66', 'A68', 'A69', 'A74', 'A80', 'A203', 'A210', 'A2403'
]

HLA_B_ALLELES = [
    'B7', 'B8', 'B13', 'B14', 'B15', 'B18', 'B27', 'B35',
    'B37', 'B38', 'B39', 'B41', 'B42', 'B44', 'B45', 'B46',
    'B47', 'B48', 'B49', 'B50', 'B51', 'B52', 'B53', 'B54',
    'B55', 'B56', 'B57', 'B58', 'B59', 'B60', 'B61', 'B62',
    'B63', 'B64', 'B65', 'B67', 'B70', 'B71', 'B72', 'B73',
    'B75', 'B76', 'B77', 'B78', 'B81', 'B82', 'B703', 'B1513',
    'B1516', 'B1517', 'B3901', 'B3902'
]

HLA_DR_ALLELES = [
    'DR1', 'DR3', 'DR4', 'DR6', 'DR7', 'DR8', 'DR9', 'DR10',
    'DR11', 'DR12', 'DR13', 'DR14', 'DR15', 'DR16', 'DR17',
    'DR18', 'DR1403', 'DR1404'
]

# assume equal distribution of HLA alleles for simplicity.
# actual frequencies can be complex and population-specific

def generateHLA(n):
    results = []

    # HLA distributions are very complex. for simplicity, we assume equal distribution of alleles. 
    for i in range(n):
        hla_a  = random.sample(HLA_A_ALLELES, 2)
        hla_b  = random.sample(HLA_B_ALLELES, 2)
        hla_dr = random.sample(HLA_DR_ALLELES, 2)
        results.append((hla_a, hla_b, hla_dr))
    return results

def generatePRA(n):
    pra_values = []
    for i in range(n):
        r = random.random()
        if r < 0.7:
            pra = round(random.betavariate(1, 8) * 100)   # mostly low PRA
        elif r < 0.9:
            pra = round(random.betavariate(2, 2) * 100)   # moderate PRA
        else:
            pra = round(random.betavariate(8, 1.5) * 100) # highly sensitized tail
        pra_values.append(pra)
    return pra_values

def assignAttributes(n):
    blood_types = generateBloodType(n)
    rh_factors  = generateRhFactor(n)
    hla_types   = generateHLA(n)
    pra_values  = generatePRA(n)

    attributes = []
    for i in range(n):
        attributes.append({
            'blood_type': blood_types[i],
            'rh_factor':  rh_factors[i],
            'hla_a':      hla_types[i][0],
            'hla_b':      hla_types[i][1],
            'hla_dr':     hla_types[i][2],
            'pra':        pra_values[i]
        })
    return attributes
