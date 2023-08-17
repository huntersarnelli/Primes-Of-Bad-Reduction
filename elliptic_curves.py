import csv

#this is a parralized algorithm from above. It works and is faster than above. 

def compute_curve_invariants(args):
    A, B = args
    E = EllipticCurve([A, B])
    
    j = E.j_invariant()
    #rank = E.rank()
    torsion_order = E.torsion_order()
    conductor = E.conductor()
    prime_factors_conductor = factor(conductor)
    
    
    segregated_primes = {}
    for prime, exponent in prime_factors_conductor:
        digit_count = len(str(prime))
        if digit_count not in segregated_primes:
            segregated_primes[digit_count] = []
        segregated_primes[digit_count].append((prime, exponent))
        
        
        
    num_p_f = len(prime_factors_conductor)
    delta = E.discriminant()
    equation = E.ainvs()
    eqn_str = "y^2 = x^3 + {}*x + {}".format(equation[3], equation[4])
    
    L = E.period_lattice()
    omega_0, omega_1 = L.basis()

    return E, j, rank, torsion_order, conductor, delta, omega_0, omega_1,eqn_str, A, B, prime_factors_conductor,num_p_f,segregated_primes



def generate_elliptic_curves(a_range, b_range, num_curves):
    curves = []
    j_invarient = []
    #e_rank = []
    tor_order = []
    conductors = []
    discrim = []
    omega0L = []
    omega1L = []
    eqn_strs = []
    AL =[]
    BL = []
    bad_red = []
    number = []
    seg_primes = []
    
    tasks = []
    while len(tasks) < num_curves:
        A = int(randint(a_range[0], a_range[1]))
        B = int(randint(b_range[0], b_range[1]))
        discriminant = -16*(4*A^3 + 27*B^2)
        if discriminant != 0:
            tasks.append((A, B))
    
    # Use as many processes as there are CPU cores
    p = Pool(cpu_count())
    results = p.map(compute_curve_invariants, tasks)
    p.close()
    p.join()
        
    for res in results:
        curves.append(res[0])
        j_invarient.append(res[1])
        #e_rank.append(res[2])
        tor_order.append(res[3])
        conductors.append(res[4])
        discrim.append(res[5])
        omega0L.append(res[6])
        omega1L.append(res[7])
        eqn_strs.append(res[8])
        AL.append(res[9])
        BL.append(res[10])
        bad_red.append(res[11])
        number.append(res[12])
        seg_primes.append(res[13])
    
    combined = list(zip(eqn_strs,j_invarient,tor_order,conductors,discrim,omega0L,omega1L,AL,BL,bad_red,number,seg_primes))
    
    

    with open('Eliptic_curves.csv', 'wb') as file:
        writer = csv.writer(file)
        writer.writerow(["eqn_strs", "j_invarient", "tor_order","conductors","discrim","omega0L","omega1L","A Coefficient","B Coefficient","bad_red","number of bad primes","seg_primes"])  # Writing headers
        writer.writerows(combined)
        
    with open('Eliptic_curves.csv', 'rb') as file:  # 'rb' for reading in binary mode in Python 2
        content = file.read()
    #print(content)


    return content,curves, j_invarient, tor_order,conductors,discrim,omega0L,omega1L,eqn_strs,A,B,bad_red,number


# Example of usage
a_range = (-50, 50)
b_range = (-50, 50)
num_curves = 1000
content, curves, j_invarient, tor_order,conductors,discrim,omega0L,omega1L,eqn_strs,A,B,bad_red,number_seg_primes = generate_elliptic_curves(a_range, b_range, num_curves)

print(content)
