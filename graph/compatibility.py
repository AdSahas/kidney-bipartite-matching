BLOOD_COMPATIBLE = {
    'O' : ['O', 'A', 'B', 'AB'],
    'A' : ['A', 'AB'],
    'B' : ['B', 'AB'],
    'AB' : ['AB']
}

def bloodCompatibility(donor_blood, recipient_blood):
    return recipient_blood in BLOOD_COMPATIBLE[donor_blood]

def rhCompatibility(donor_rh, recipient_rh):
    # Rh+ can only donate to Rh+.
    if donor_rh == '+':
        return recipient_rh == '+'
    
    else:
        return True # Rh- is a universal Rh factor donor. 
    

def HLAMatchScore(donor, recipient):
    score = 0
    for d_alleles, r_alleles in zip(
        [donor['hla_a'],     donor['hla_b'],     donor['hla_dr']],
        [recipient['hla_a'], recipient['hla_b'], recipient['hla_dr']]
    ):
        score += len(set(d_alleles) & set(r_alleles))
    return score  # [0, 6] where 6 is a perfect match

def crossmatch_compatible(donor, recipient):
    score = HLAMatchScore(donor, recipient)
    pra = recipient['pra']
    rejection_chance = (pra / 100) * (1 - score / 6)
    return rejection_chance 

def kidneyCompatibility(donor, recipient):
    return (
        bloodCompatibility(donor['blood_type'], recipient['blood_type']) and
        rhCompatibility(donor['rh_factor'], recipient['rh_factor'])
    )