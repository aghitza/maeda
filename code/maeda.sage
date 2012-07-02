def doit(weights):
    """
    """
    ni = os.nice(19)
    for irred in maeda_parallel(weights):
        print irred

def doit_consec(weights):
    """
    """
    ni = os.nice(19)
    for irred in maeda_parallel_consec(weights):
        print irred
              
@parallel(ncpus=6)
def maeda_parallel(k):
    stk = str(k)
    filename = 'data/' + '0'*(5-len(stk)) + stk
    lockfilename = filename + '.lock'
    if os.path.exists(filename) or os.path.exists(lockfilename):
        irred = None
    else:
        import time as time
        import socket as socket
        lockfile = open(lockfilename, 'w')
        lockfile.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))
        lockfile.write(' ' + socket.gethostname())
        lockfile.write('\n')
        lockfile.close()
        irred = maeda_modular(weight=k, verbose=True, filename=filename)
        os.remove(lockfilename)
    return irred

@parallel(ncpus=6)
def maeda_parallel_consec(k):
    stk = str(k)
    filename = 'data-consec/' + '0'*(5-len(stk)) + stk
    lockfilename = filename + '.lock'
    if os.path.exists(filename) or os.path.exists(lockfilename):
        irred = None
    else:
        import time as time
        import socket as socket
        lockfile = open(lockfilename, 'w')
        lockfile.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))
        lockfile.write(' ' + socket.gethostname())
        lockfile.write('\n')
        lockfile.close()
        irred = maeda_modular_consec(weight=k, verbose=True, filename=filename)
        os.remove(lockfilename)
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

def maeda_modular(weight, PRIME_BOUND=2^17, verbose=True, filename=None):
    import time as time
    time0 = time.time()
    count = 0
    dim = dimension_cusp_forms(1, weight)
    prec = 2*(dim + 2)
    #b = CuspForms(1, weight).q_expansion_basis(prec=prec)
    b = victor_miller_basis(weight, prec=prec, cusp_only=True)
    time1 = time.time()
    time_str = ''
    st = "# basis computed in              %9.3f seconds"%(time1-time0)
    print(st)
    time_str += st + '\n'
    M = hecke_operator_on_basis(b, 2, weight)
    time2 = time.time()
    st = "# matrix of T_2 computed in      %9.3f seconds"%(time2-time1)
    print(st)
    time_str += st + '\n'
    type1_prime = None
    type2_prime = None
    type3_prime = None
    while ((type1_prime is None) or (type2_prime is None) or
        (type3_prime is None)):
        count += 1
        p = random_prime(PRIME_BOUND)
        K = GF(p)
        Mp = M.change_ring(K)
        f = Mp.charpoly()
        if (type1_prime is None) and f.is_irreducible():
            type1_prime = p
            type1_poly = f
            st = "# type 1 found in                %9.3f seconds, after %4d tries"%(time.time()-time2, count)
            print(st)
            time_str += st + '\n'
            if (type3_prime is None) and is_prime(f.degree()):
                type3_prime = p
                type3_poly = f
                st = "# type 3 found in                %9.3f seconds, after %4d tries"%(time.time()-time2, count)
                print(st)
                time_str += st + '\n'
                continue
        if not f.is_squarefree():
            continue
        fact = f.factor()
        lst = sorted([g[0].degree() for g in fact])
        if (type2_prime is None) and is_type_II(lst):
            type2_prime = p
            type2_poly = f
            st = "# type 2 found in                %9.3f seconds, after %4d tries"%(time.time()-time2, count)
            print(st)
            time_str += st + '\n'
        if (type3_prime is None) and is_type_III(lst):
            type3_prime = p
            type3_poly = f
            st = "# type 3 found in                %9.3f seconds, after %4d tries"%(time.time()-time2, count)
            print(st)
            time_str += st + '\n'
    if not filename is None:
        file = open(filename, 'w')
        file.write("# weight %s\n" %weight)
        file.write("# total time taken (in seconds): %9.3f\n" %(time.time() - time0))
        file.write(time_str)
        file.write("# %s\n" %version())
        file.write("# %s\n" %(' '.join(os.uname()[:3]) + ' ' + os.uname()[4]))
        file.write("# %s\n" %(time.asctime()))
        file.write("# type 1 prime and polynomial\n")
        file.write("%s\n" %type1_prime)
        file.write("%s\n" %type1_poly)
        file.write("# type 2 prime and polynomial\n")
        file.write("%s\n" %type2_prime)
        file.write("%s\n" %type2_poly)
        file.write("# type 3 prime and polynomial\n")
        file.write("%s\n" %type3_prime)
        file.write("%s\n" %type3_poly)
        file.close()     
    return True

def maeda_modular_consec(weight, PRIME_BOUND=2^20, verbose=True, filename=None):
    import time as time
    time0 = time.time()
    count = 0
    dim = dimension_cusp_forms(1, weight)
    prec = 2*(dim + 2)
    #b = CuspForms(1, weight).q_expansion_basis(prec=prec)
    b = victor_miller_basis(weight, prec=prec, cusp_only=True)
    time1 = time.time()
    time_str = ''
    st = "# basis computed in              %9.3f seconds"%(time1-time0)
    print(st)
    time_str += st + '\n'
    M = hecke_operator_on_basis(b, 2, weight)
    time2 = time.time()
    st = "# matrix of T_2 computed in      %9.3f seconds"%(time2-time1)
    print(st)
    time_str += st + '\n'
    type1_prime = None
    type2_prime = None
    type3_prime = None
    p = 1
    while ((type1_prime is None) or (type2_prime is None) or
        (type3_prime is None)):
        count += 1
        p = next_prime(p)
        K = GF(p)
        Mp = M.change_ring(K)
        f = Mp.charpoly()
        if (type1_prime is None) and f.is_irreducible():
            type1_prime = p
            type1_poly = f
            st = "# type 1 found in                %9.3f seconds, after %4d tries"%(time.time()-time2, count)
            print(st)
            time_str += st + '\n'
            if (type3_prime is None) and is_prime(f.degree()):
                type3_prime = p
                type3_poly = f
                st = "# type 3 found in                %9.3f seconds, after %4d tries"%(time.time()-time2, count)
                print(st)
                time_str += st + '\n'
                continue
        if not f.is_squarefree():
            continue
        fact = f.factor()
        lst = sorted([g[0].degree() for g in fact])
        if (type2_prime is None) and is_type_II(lst):
            type2_prime = p
            type2_poly = f
            st = "# type 2 found in                %9.3f seconds, after %4d tries"%(time.time()-time2, count)
            print(st)
            time_str += st + '\n'
        if (type3_prime is None) and is_type_III(lst):
            type3_prime = p
            type3_poly = f
            st = "# type 3 found in                %9.3f seconds, after %4d tries"%(time.time()-time2, count)
            print(st)
            time_str += st + '\n'
    if not filename is None:
        file = open(filename, 'w')
        file.write("# weight %s\n" %weight)
        file.write("# total time taken (in seconds): %9.3f\n" %(time.time() - time0))
        file.write(time_str)
        file.write("# %s\n" %version())
        file.write("# %s\n" %(' '.join(os.uname()[:3]) + ' ' + os.uname()[4]))
        file.write("# %s\n" %(time.asctime()))
        file.write("# type 1 prime and polynomial\n")
        file.write("%s\n" %type1_prime)
        file.write("%s\n" %type1_poly)
        file.write("# type 2 prime and polynomial\n")
        file.write("%s\n" %type2_prime)
        file.write("%s\n" %type2_poly)
        file.write("# type 3 prime and polynomial\n")
        file.write("%s\n" %type3_prime)
        file.write("%s\n" %type3_poly)
        file.close()     
    return True

def is_type_II(lst):
    mylst = copy(lst)
    if mylst.count(2) != 1:
        return False
    mylst.remove(2)
    for l in mylst:
        if (l % 2) == 0:
            return False
    return True

def is_type_III(lst):
    return (is_prime(lst[-1]) and (lst[-1] > sum(lst)/2))

def check_results(pathname):
    for fname in sorted(os.listdir(pathname)):
        print("%s"%fname)
        fullname = pathname + fname
        resfile = open(fullname, 'r')
        for _ in range(11):
            resfile.readline()
        p = ZZ(resfile.readline().strip('\n'))
        R.<x> = GF(p)[]
        poly = resfile.readline().strip('\n')
        f = R(poly)
        if not f.is_irreducible():
            print("%s: type 1 is reducible"%fname)
        resfile.readline()
        p = ZZ(resfile.readline().strip('\n'))
        R.<x> = GF(p)[]
        poly = resfile.readline().strip('\n')
        f = R(poly)
        if not f.is_squarefree():
            print("%s: type 2 is not squarefree"%fname)
        else:
            fact = f.factor()
            lst = sorted([g[0].degree() for g in fact])
            if not is_type_II(lst):
                print("%s: type 2 is not the right shape"%fname)
        resfile.readline()
        p = ZZ(resfile.readline().strip('\n'))
        R.<x> = GF(p)[]
        poly = resfile.readline().strip('\n')
        f = R(poly)
        if not f.is_squarefree():
            print("%s: type 3 is not squarefree"%fname)
        else:
            fact = f.factor()
            lst = sorted([g[0].degree() for g in fact])
            if not is_type_III(lst):
                print("%s: type 3 is not the right shape"%fname)
        resfile.close()

def get_stats(pathname):
    type1 = dict()
    type2 = dict()
    type3 = dict()
    for fname in sorted(os.listdir(pathname)):
        print("%s"%fname)
        if ".lock" in fname:
            continue
        k = ZZ(fname.lstrip('0'))
        d = dimension_cusp_forms(1, k)
        fullname = pathname + fname
        resfile = open(fullname, 'r')
        for _ in range(4):
            resfile.readline()
        for _ in range(3):
            primeline = resfile.readline().strip('\n')
            strlist = primeline.split(' ')
            tp = ZZ(strlist[2])
            tries = ZZ(strlist[-2])
            if tp == 1:
                type1[k] = tries/d*1.0
            elif tp == 2:
                type2[k] = tries
            elif tp == 3:
                type3[k] = tries
            else:
                raise RuntimeError("booom")
        resfile.close()
    return type1, type2, type3

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
