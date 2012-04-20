def doit(weights):
    """
    """
    ni = os.nice(19)
    for irred in maeda_parallel(weights):
        print irred
              
@parallel(ncpus=10)
def maeda_parallel(k):
    stk = str(k)
    filename = 'data/' + '0'*(5-len(stk)) + stk
    if os.path.exists(filename):
        irred = None
    else:
        irred = maeda(weight=k, verbose=True, filename=filename)
    return irred

def maeda(weight, verbose=True, filename=None):
    """

    INPUT:

    - ``weight`` -- integer
    - ``verbose`` (default: True) -- indicates whether to give some
      information about what is happening
    - ``filename`` (default: None) -- if not None, a string giving the file
      name under which to save the results
                            
    EXAMPLES::

    """
    import time as time
    #mem0 = get_memory_usage()
    dim = dimension_cusp_forms(1, weight)
    prec = 2 * (dim + 2)
    if verbose:
        print "%s: computing the basis of the space..." %weight
    time0 = time.time()
    b = CuspForms(1, weight).q_expansion_basis(prec=prec)
    time1 = time.time()
    if verbose:
        print "%s: ... done in %s seconds" %(weight, (time1 - time0))
        print "%s: computing the Hecke matrix..." %weight
    M = hecke_operator_on_basis(b, 2, weight)
    #mem1 = get_memory_usage()
    del b
    time2 = time.time()
    if verbose:
        print "%s: ... done in %s seconds" %(weight, (time2 - time1))
        print "%s: computing the characteristic polynomial..." %weight
    f = M.charpoly()
    #mem2 = get_memory_usage()
    del M
    time3 = time.time()
    if verbose:
	    print "%s: ... done in %s seconds" %(weight, (time3 - time2))
	    print "%s: testing for irreducibility..." %weight
    isirred = f.is_irreducible()
    #mem3 = get_memory_usage()
    time4 = time.time()
    if verbose:
        print "%s: ... done in %s seconds" %(weight, (time4 - time3))
        #print "%s: maximum memory usage: %s" %(weight, max([mem3, mem2, mem1]) - mem0)
    if not filename is None:
        file = open(filename, 'w')
        file.write("%s\n" %f)
        del f
        file.write("%s\n" %isirred)
        file.write("# weight %s\n" %weight)
        file.write("# total time taken (in seconds): %s\n" %(time4 - time0))
        file.write("#     find basis:                %s\n" %(time1 - time0))
        file.write("#     compute matrix:            %s\n" %(time2 - time1))
        file.write("#     find charpoly:             %s\n" %(time3 - time2))
        file.write("#     irred test:                %s\n" %(time4 - time3))
        #file.write("# maximum memory usage (in Mb): %s\n" %(max([mem3, mem2, mem1]) - mem0))
        #file.write("#     basis and matrix: %s\n" %(mem1 - mem0))
        #file.write("#     matrix and poly:  %s\n" %(mem2 - mem0))
        #file.write("#     poly:             %s\n" %(mem3 - mem0))
        file.write("# %s\n" %version())
        file.write("# %s\n" %(' '.join(os.uname()[:3]) + ' ' + os.uname()[4]))
        file.write("# %s\n" %(time.asctime()))
        file.close()     
    else:
        del f
    return isirred         

def maeda_modular(weight, basis, verbose=True, filename=None):
    import time as time
    count = 0
    while True:
        count += 1
        p = random_prime(2^20)
        K = GF(p)
        bp = [f.change_ring(K) for f in b]
        M = hecke_operator_on_basis(bp, 2, weight)
        del bp
        f = M.charpoly()
        del M
        isirred = f.is_irreducible()
        del f
        if isirred:
            return p, count

def maeda_modp(weight, p, verbose=True, filename=None):
    import time as time
    dim = dimension_cusp_forms(1, weight)
    prec = 2 * (dim + 2)
    if verbose:
        print "%s: computing the basis of the space..." %weight
    time0 = time.time()
    b = CuspForms(1, weight).q_expansion_basis(prec=prec)
    b = [f.change_ring(GF(p)) for f in b]
    time1 = time.time()
    if verbose:
        print "%s: ... done in %s seconds" %(weight, (time1 - time0))
        print "%s: computing the Hecke matrix..." %weight
    M = hecke_operator_on_basis(b, 2, weight)
    #mem1 = get_memory_usage()
    del b
    time2 = time.time()
    if verbose:
        print "%s: ... done in %s seconds" %(weight, (time2 - time1))
        print "%s: computing the characteristic polynomial..." %weight
    f = M.charpoly()
    #mem2 = get_memory_usage()
    del M
    time3 = time.time()
    if verbose:
	    print "%s: ... done in %s seconds" %(weight, (time3 - time2))
	    print "%s: testing for irreducibility..." %weight
    isirred = f.is_irreducible()
    #mem3 = get_memory_usage()
    time4 = time.time()
    if verbose:
        print "%s: ... done in %s seconds" %(weight, (time4 - time3))
        #print "%s: maximum memory usage: %s" %(weight, max([mem3, mem2, mem1]) - mem0)
    return isirred         
