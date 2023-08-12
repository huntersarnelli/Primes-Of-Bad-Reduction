
# coding: utf-8

# In[55]:


def buchberger_with_integers(I):
    """
    Modified Buchberger's algorithm to compute a Gröbner basis for ideal I 
    and to capture all integer coefficients that appear during the computation.
    """
    def capture_coefficients(poly, captured):
        """Utility function to extract and store coefficients of a polynomial."""
        for coeff in poly.coefficients():
            if coeff.parent() in [QQ, ZZ]:
                captured.add(Integer(round(coeff)))
            else:
                captured.add(Integer(round(coeff.real())))
                captured.add(Integer(round(coeff.imag())))

    basis = list(I.gens())
    pairs = [(basis[i], basis[j]) for i in range(len(basis)) for j in range(i+1, len(basis))]
    
    captured_integers = set()

    for poly in basis:
        capture_coefficients(poly, captured_integers)

    while pairs:
        f, g = pairs.pop()

        capture_coefficients(f, captured_integers)
        capture_coefficients(g, captured_integers)
        
        S = I.reduce(f.lc() * g - g.lc() * f)
        
        capture_coefficients(S, captured_integers)

        if not S.is_zero():
            for h in basis:
                pairs.append((S, h))
            
            basis.append(S)
            pairs = list(set(pairs))

    for poly in basis:
        capture_coefficients(poly, captured_integers)

    return basis, captured_integers


def ideal_intersection(I, J):
    """
    Compute the intersection of two ideals I and J using Gröbner bases.
    """
    def capture_coefficients_from_ideal(ideal, captured):
        """Utility function to extract and store coefficients from an ideal."""
        for poly in ideal.gens():
            for coeff in poly.coefficients():
                if coeff.parent() in [QQ, ZZ]:
                    captured.add(Integer(coeff))
                else:
                    captured.add(Integer(round(coeff.real())))
                    captured.add(Integer(round(coeff.imag())))

    if I.ring() != J.ring():
        raise ValueError("The ideals must be in the same polynomial ring")

    R = I.ring()
    variables = ['t'] + list(R.variable_names())
    S = PolynomialRing(R.base_ring(), variables, order='lex')
    t = S.gen(0)

    I_gens_in_S = [S(f) for f in I.gens()]
    J_gens_in_S = [S(f) for f in J.gens()]
    K_gens = [t*f for f in I_gens_in_S] + [(1-t)*g for g in J_gens_in_S]
    K = S.ideal(K_gens)
    
    print('here is I_gens', I_gens_in_S)
    print('here is J_gens', J_gens_in_S)

    # Capture all intermediary calculations
    G, integers = buchberger_with_integers(K)
    print(integers)
    G1 = K.groebner_basis()

    intersection_gens = [f for f in G1 if t not in f.variables()]
    
    intersection_ideal = R.ideal(intersection_gens)

    # Capture coefficients from the final intersection ideal
    capture_coefficients_from_ideal(intersection_ideal, integers)

    return intersection_ideal, integers



def ideal_quotient(I, J):
    """
    Compute the quotient of two ideals I and J.
    """
    if I.ring() != J.ring():
        raise ValueError("The ideals must be in the same polynomial ring")

    R = I.ring()
    quotient_ideals = []
    total_integers = set()

    for g in J.gens():
        # Capture intersection and all intermediary calculations
        intersection, ints = ideal_intersection(I, R.ideal(g))
        total_integers.update(ints)
        
        scaled_ideal = R.ideal([(1/g)*f for f in intersection.gens()])
        quotient_ideals.append(scaled_ideal)

    if len(quotient_ideals) == 1:
        return quotient_ideals[0], total_integers

    final_ideal, ints = ideal_intersection(quotient_ideals[0], quotient_ideals[1])
    total_integers.update(ints)

    for i in range(2, len(quotient_ideals)):
        final_ideal, ints = ideal_intersection(final_ideal, quotient_ideals[i])
        total_integers.update(ints)

    return final_ideal, total_integers

R.<x,y,z> = ZZ[]
f = f = x^2*y^2 - y^4 + x^3*y + x^4 +z^4
fx = f.derivative(x)
fy = f.derivative(y)
fz = f.derivative(z)
I = R.ideal([f,fx,fy,fz])
J = R.ideal(x^30,y^30,z^30)

G_custom, integers1 = buchberger_with_integers(I)

#print("Custom Gröbner Basis:", G_custom)
print("Integers captured:", integers1)

intersection2, integers2 = ideal_intersection(I, J)

#print("Intersection:", intersection2)
print("Integers captured:", integers2)

quotient, integers3 = ideal_quotient(I, J)

#print("Ideal quotient:", quotient)
print("Integers captured:", integers3)

#11342685 = 3 ⋅ 5 ⋅ 756179

