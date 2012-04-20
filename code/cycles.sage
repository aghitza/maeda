def is_type_II(g):
    lst = sorted([len(c) for c in g.cycle_tuples()])
    if lst.count(2) != 1:
        return False
    lst.remove(2)
    for l in lst:
        if (l % 2) == 0:
            return False
    return True

def is_type_III(g, d):
    lst = uniq([len(c) for c in g.cycle_tuples()])
    for l in lst:
        if is_prime(l) and (l > d/2):
            return True
    return False

def proportion_of_type_II(d, method="other"):
    if method == "symmetric":
        S = SymmetricGroup(d)
        count = 0
        for g in S:
            if is_type_II(g):
                count += 1
        return count/factorial(d)
    if method == "permutation":
        S = Permutations(d)
        count = 0
        for g in S:
            if is_type_II(g):
                count += 1
        return count/factorial(d)
    P = Partitions(d-2, parts_in=[1,3..(d-2)])
    sum = 0
    for p in P:
        l = p.to_list()
        fact = 1
        for j in range(min(l), max(l)+1):
            fact *= factorial(l.count(j))
        sum += 1/(prod(l)*fact)
    return sum/2

def proportion_of_type_III(d, method="other"):
    if method == "symmetric":
        S = SymmetricGroup(d)
        count = 0
        for g in S:
            if is_type_III(g, d):
                count += 1
        return count/factorial(d)
    if method == "permutation":
        S = Permutations(d)
        count = 0
        for g in S:
            if is_type_III(g, d):
                count += 1
        return count/factorial(d)
    sum2 = 0
    for q in prime_range((d+1)/2, d+1):
        P = Partitions(d-q)
        sum1 = 0
        for p in P:
            l = p.to_list() + [q]
            fact = 1
            for j in range(min(l), max(l)+1):
                fact *= factorial(l.count(j))
            sum1 += 1/(prod(l)*fact)
        sum2 += sum1
    return sum2
