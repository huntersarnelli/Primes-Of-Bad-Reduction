import csv

#this is a parralized algorithm from above. It works and is faster than above. 

def compute_curve_invariants(args):
    A, B = args
    E = EllipticCurve([A, B])
    
    j = E.j_invariant()
    #rank = E.rank()
    torsion_order = E.torsion_order()
    conductor = E.conductor()
    delta = E.discriminant()
    equation = E.ainvs()
    eqn_str = "y^2 = x^3 + {}*x + {}".format(equation[3], equation[4])
    
    L = E.period_lattice()
    omega_0, omega_1 = L.basis()

    return E, j, rank, torsion_order, conductor, delta, omega_0, omega_1,eqn_str, A, B



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
    
    combined = list(zip(eqn_strs,j_invarient,tor_order,conductors,discrim,omega0L,omega1L,AL,BL))
    
    

    with open('Eliptic_curves.csv', 'wb') as file:
        writer = csv.writer(file)
        writer.writerow(["eqn_strs", "j_invarient", "tor_order","conductors","discrim","omega0L","omega1L","A Coefficient","B Coefficient"])  # Writing headers
        writer.writerows(combined)
        
    with open('Eliptic_curves.csv', 'rb') as file:  # 'rb' for reading in binary mode in Python 2
        content = file.read()
    #print(content)


    return content,curves, j_invarient, tor_order,conductors,discrim,omega0L,omega1L,eqn_strs,A,B


# Example of usage
a_range = (0, 50)
b_range = (0, 50)
num_curves = 100
content, curves, j_invarient, tor_order,conductors,discrim,omega0L,omega1L,eqn_strs,A,B = generate_elliptic_curves(a_range, b_range, num_curves)

print(content)
