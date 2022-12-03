#!/usr/bin/env sage
from itertools import combinations_with_replacement, chain, combinations, permutations
import math
from itertools import product as gen_product

# WARNING both n and t are SYMBOLIC variables!!!
var('n t')

def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(2, len(s)+1))

def potential_homomorphism(temp1, temp2):
    blocks1 = sorted(temp1.block_sizes())
    blocks2 = sorted(temp2.block_sizes())

    # First cannot have more blocks than second
    if len(blocks1) > len(blocks2):
        return False

    if blocks1 == blocks2:
        return False

    return True

    # Additionally the blocks of the first template must be smaller (non-strictly)
    if any(b1 > b2 for b1, b2 in zip(blocks1, blocks2)):
        return False

    # And at least one smaller
    if all(b1 == b2 for b1, b2 in zip(blocks1, blocks2)):
        return False

    return True

def generate_mapping(images, used, partial=[]):
    if len(images) == 0:
        yield partial
    else:
        for y in images[0]:
            if y in used:
                continue
            yield from generate_mapping(images[1:], used | {y}, partial=partial+[y])
        
def permutations_with_repetition(iterable, length, used=set(), partial=[]):
    if length == 0:
        yield partial
    unused = set(iterable)-used
    if len(unused) == length:
        for perm in permutations(unused):
            yield partial + list(perm)
    elif len(unused) < length:
        for elm in iterable:
            yield from permutations_with_repetition(iterable, length-1, used|{elm}, partial+[elm])

def count_homomorphisms(temp1, temp2):
    homs = set()
    blocks1 = temp1.blocks()
    blocks2 = temp2.blocks()

    for mapping in permutations_with_repetition(list(range(len(blocks1))), len(blocks2)):
        mapping = [blocks1[i] for i in mapping]
        # Compute potential images of every point in temp1
        images = [set(temp2.ground_set()) for _ in range(len(temp1.blocks()))]

        # And restrict them to the respective parts
        for b1, b2 in zip(mapping, blocks2):
            for x in b1:
                images[x] &= set(b2)
                
        # Try all valid mappings
        for func in generate_mapping(images,set()):
            # Each non-isomorphic homomorphism is defined by the set of unused points,
            # together with the blocks of mod_
            homs |= {tuple(sorted(tuple(sorted(func[x] for x in block)) for block in blocks1))}
    #print(homs)
    return len(homs)

def comb_number(i,j):
    return product(n-k for k in range(i,j))

def generate_templates(dim):
    templates = []
    for num in range(dim+1):
        for template in hypergraphs.nauty(num, dim*dim, multiple_sets=False,
                                          set_max_size=dim,set_min_size=1):
            structure = IncidenceStructure(template)
            valid = True
            for subset in powerset(template):
                intersection = list(x for x in structure.ground_set() if all(x in s for s in subset))
                if len(intersection) >= sum(len(s) for s in subset) - (len(subset) - 1)*(dim+1):
                    valid = False
                    break

            if valid:
                templates.append(structure)
    return templates
                

def compute_polynomial(dim):
    # Compute all templates
    templates = generate_templates(dim)

    dim_templates = [[] for _ in range(dim + 1)]
    
    for i, temp in enumerate(templates):
        temp_dim = dim - sum((dim+1) - bl for bl in temp.block_sizes())
        dim_templates[temp_dim].append(i)

    #for d in range(dim+1):
    #    print("Dimenze {}:".format(d))
    #    for i in dim_templates[d]:
    #        print(templates[i].blocks())

    
    contributions = [0]*len(templates)
    polynomial = 0

    for d in reversed(range(dim + 1)):
        for i in dim_templates[d]:
            template = templates[i]
            # Basic contribution
            cont = (-1)^(len(template.blocks()))

            for j in chain.from_iterable(dim_templates[d2] for d2 in range(d+1, dim+1)):
                other = templates[j]

                # Consider only size matches
                if not potential_homomorphism(template, other):
                    continue

                # subtract the contribution of other multiplied by the number of homomorphisms
                cont -= (comb_number(len(template.ground_set()),len(other.ground_set()))*
                        count_homomorphisms(template, other)*
                        template.automorphism_group().cardinality()/
                        other.automorphism_group().cardinality()*contributions[j])
                
            #print(i, template.blocks(), cont)

            # Save the actual contribution
            contributions[i] = cont

            # And update polynomial
            polynomial += (cont*comb_number(0,len(template.ground_set()))/
                    template.automorphism_group().cardinality()*t^d)

    # Profit
    return polynomial


