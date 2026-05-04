from attributes import assignAttributes
from compatibility import kidneyCompatibility
from compatibility import crossmatch_compatible, HLAMatchScore, kidneyCompatibility, bloodCompatibility, rhCompatibility

INF = 1e9

# create a bipartite graph.
def build_graph(n_donors, n_recipients):
    donors = assignAttributes(n_donors)
    recipients = assignAttributes(n_recipients)
    return donors, recipients

# greedy matching implementation
def GreedyMatching(donors, recipients, rej_chance):
    matches = [-1] * len(recipients)
    for i, donor in enumerate(donors):
        for r in range(len(recipients)):
            if matches[r] == -1 and kidneyCompatibility(donor, recipients[r]) and crossmatch_compatible(donor, recipients[r]) < rej_chance:
                matches[r] = i
                break
    return matches




# implement kuhn's algorithm
def kuhnDFS(donorIdx, donors, recipients, matches, visited, rej_chance):

    # keep track of donors
    donor = donors[donorIdx]

    # for every recipient
    for r in range(len(recipients)):

        # check if they are visited. if so, skip.
        if visited[r]:
            continue

        # otherwise, check compatibility. 
        if kidneyCompatibility(donor, recipients[r]) and crossmatch_compatible(donor, recipients[r]) < rej_chance:

            # check whether they are matched or not.
            visited[r] = True
            if matches[r] == -1:

                # match if not matched
                matches[r] = donorIdx  # store index
                return True
            
            # if matched, can an alternative be found? use recursive DFS. 
            if kuhnDFS(matches[r], donors, recipients, matches, visited, rej_chance):
                matches[r] = donorIdx  # store index
                return True
            
    # if no match found, return false. 
    return False

def kuhn(donors, recipients, rej_chance):

    # keep track of which recipients are matched. -1 if unmatched. 
    matches = [-1] * len(recipients)

    # for every donor, try to find a match with DFS. 
    for i in range(len(donors)):
        visited = [False] * len(recipients)
        kuhnDFS(i, donors, recipients, matches, visited, rej_chance)
    return matches



# Hungarian implementation - jonker volgenant verison for better performance.
def getCost(donor, recipient, rej_chance):

    # only calculate costs for blood and Rh comparible pairs.
    if not bloodCompatibility(donor['blood_type'], recipient['blood_type']):
        return INF
    if not rhCompatibility(donor['rh_factor'], recipient['rh_factor']):
        return INF
    
    # using rejection chance and hla score to calculate cost. 
    hla_score = HLAMatchScore(donor, recipient)
    rejection_chance = crossmatch_compatible(donor, recipient)

    # if the rejection chance is above threshold, set it to inf - cannot be considered.
    if rejection_chance >= rej_chance:
        return INF
    return rejection_chance * (1 - hla_score / 6)

def dijkstra(donor, donors, recipients, u, v, matchD, matchR, rej_chance):
    n = len(recipients)
    dist = [INF] * n
    visited = [False] * n

    # keep track of previous donor and recipient in the path
    prevD = [-1] * n
    prevR = [-1] * n

    # initialize distances from starting donor
    for r in range(n):
        cost = getCost(donors[donor], recipients[r], rej_chance)
        if cost < INF:
            dist[r] = cost - u[donor] - v[r]

    while True:
        # pick unvisited recipient with smallest dist
        rMin = -1

        for r in range(n):
            if not visited[r]:
                if rMin == -1 or dist[r] < dist[rMin]:
                    rMin = r

        # no augmenting path exists
        if rMin == -1 or dist[rMin] == INF:
            return None
        
        visited[rMin] = True
        # if recipient is unmatched, we found an augmenting path
        if matchR[rMin] == -1:
            return rMin, dist, prevD, prevR
        
        # else, we need to relax edges from the matched donor
        dNext = matchR[rMin]
        for r in range(n):
            if not visited[r]:
                cost = getCost(donors[dNext], recipients[r], rej_chance)
                if cost < INF:
                    newDist =  dist[rMin] + cost - u[dNext] - v[r]
                    if newDist < dist[r]:
                        dist[r] = newDist
                        prevD[r] = dNext
                        prevR[r] = rMin


def jonkerVolgenant(donors, recipients, rej_chance):

    '''
    this is a very complicated and challenging algorithm, which took me several hours
    to understand. it was not intuitive in any way, and an understanding of
    linear programming is required to understand the update rules. 
    i had to use claude to explain the algorithm through examples,
    and learn about the duality involved. in the appendix of my report,
    this is explored in brief detail. 
    another resource: https://medium.com/@rajneeshtiwari_22870/linear-assignment-problem-in-metric-learning-for-computer-vision-eba7d637c5d4
    '''
    nD = len(donors)
    nR = len(recipients)
    u = [0.0] * nD   # donor potentials
    v = [0.0] * nR   # recipient potentials
    matchD = [-1] * nD # matchD[donor] = recipient
    matchR = [-1] * nR # matchR[recipient] = donor

    # for every donor, find an augmenting path to unmatched recipient.
    for donor in range(nD):
        # run dijkstra from this unmatched donor
        result = dijkstra(donor, donors, recipients, u, v, matchD, matchR, rej_chance)

        # no augmenting path found, skip this donor
        if result is None:
            continue

        rFree, dist, prevD, prevR = result

        # update potentials
        for r in range(nR):
            if dist[r] < INF:
                v[r] += dist[r] - dist[rFree]
        for d in range(nD):
            if matchD[d] != -1:
                u[d] += dist[matchD[d]] - dist[rFree]
        u[donor] += dist[rFree]

        # trace back path and augment
        r = rFree
        while r != -1:
            d = prevD[r]
            rPrev = prevR[r]
            if d == -1:
                matchR[r] = donor
                matchD[donor] = r
            else:
                matchR[r] = d
                matchD[d] = r
            r = rPrev

    return [(d, matchD[d]) for d in range(nD) if matchD[d] != -1]
