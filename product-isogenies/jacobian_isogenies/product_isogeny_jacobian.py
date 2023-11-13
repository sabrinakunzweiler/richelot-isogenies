import time

from sage.all import ZZ, inverse_mod
#from utilities.strategy import optimised_strategy
from jacobian_isogenies.formulae_jacobian import *
from jacobian_isogenies.doubling_jacobian import *
from jacobian_isogenies.gluing_isogeny_jacobian import *
from jacobian_isogenies.generic_isogeny_jacobian import *
from jacobian_isogenies.splitting_isogeny_jacobian import *



def optimised_strategy_kummer(n, mul_c=1):
    """
    Algorithm 60: https://sike.org/files/SIDH-spec.pdf
    
    Shown to be appropriate for (l,l)-chains in 
    https://ia.cr/2023/508
    
    Note: the costs we consider are:
       eval_c: the cost of one isogeny evaluation
       mul_c:  the cost of one element doubling

    Note: this is a slightly different version, to avoid doubling on the second Kummer surface.
    """

    eval_c = 1.000
    mul_c  = mul_c

    S = {1:[]}
    C = {1:0 }
    for i in range(2, n+1):
        b, cost = min(((b, C[i-b] + C[b] + b*mul_c + (i-b)*eval_c) for b in range(1,i)), key=lambda t: t[1])
        S[i] = [b] + S[i-b] + S[b]
        C[i] = cost

    return S[n]

def isogeny_chain_jacobian(kernel, n):
        """
        We compute a (2^n,2^n)-chain.
        The input is expected to be a (2^(n+1),2^(n+1))-kernel.
        """

        # Extract points from kernel
        (P1, Q1, P2, Q2) = kernel

        Tp1 = [P1, P2]
        Tp2 = [Q1, Q2]

        # Store chain of (2,2)-isogenies
        # we store special data associated to Type1-Type2 isogenies
        isogeny_chain = []
        phi = 0

        # Bookkeeping for optimal strategy
        strat_idx = 0
        level = [0]
        ker = [Tp1, Tp2]
        kernel_elements = [ker]
        strategy = optimised_strategy_kummer(n, mul_c=0.4)
        t_double_ell = 0
        t_double_jac = 0
        t_glue_points = 0
        t_jac_points = 0
        t_compute_glue = 0
        t_compute_jacobian_constants = 0
        t_compute_split = 0
        for k in range(n):
            prev = sum(level)
            ker = kernel_elements[-1]
            t0 = time.process_time()
            while prev != (n - 1 - k):
                level.append(strategy[strat_idx])

                # Perform the doublings
                Tp1 = double(ker[0],strategy[strat_idx],phi)
                Tp2 = double(ker[1],strategy[strat_idx],phi)

                ker = [Tp1, Tp2]

                # Update kernel elements and bookkeeping variables
                kernel_elements.append(ker)
                prev += strategy[strat_idx]
                strat_idx += 1

            # Compute the codomain from the 4-torsion
            if k == 0:
                t_double_ell = t_double_ell +  time.process_time() - t0
                t0 = time.process_time()
                phi = gluing_jacobian_isogeny(ker)
                t_compute_glue = t_compute_glue + time.process_time() - t0
            elif k == n-1:
                #is this always 4-torsion?
                t0 = time.process_time()
                phi = splitting_jacobian_isogeny(ker, phi)
                #isogeny_chain.append(splitting_isogeny)
                t_compute_split = t_compute_split + time.process_time() - t0
            else:
                t_double_jac = t_double_jac +  time.process_time() - t0
                t0 = time.process_time()
                phi = jacobian_isogeny(ker, phi)
                t_compute_jacobian_constants = t_compute_jacobian_constants + time.process_time() - t0
                        
            # Update the chain of isogenies
            isogeny_chain.append(phi)

            # Remove elements from list
            kernel_elements.pop()
            level.pop()

            # Push through points for the next step
            t0 = time.process_time()
            if k == 0:
                kernel_elements = [[apply_gluing(ker[0],phi), apply_gluing(ker[1],phi)] for ker in kernel_elements]
                t_glue_points = t_glue_points + time.process_time() - t0
                #delete later
                #print("new 2-torsion with type1:", apply_gluing(ker[0],phi), apply_gluing(ker[0],phi))
            elif k != n-1:
                kernel_elements = [(apply_jacobian_isogeny(ker[0],phi), apply_jacobian_isogeny(ker[1],phi)) for ker in kernel_elements]
                t_jac_points = t_jac_points + time.process_time() - t0
        """ some timings
        print("doublings on elliptic curve:", t_double_ell)
        print("doublings on the jacobian:", t_double_jac)
        print("compute the gluing constants:", t_compute_glue)
        print("compute the jacobian constants:", t_compute_jacobian_constants)
        print("compute the splitting constants:", t_compute_split)
        print("pushing points through gluing:", t_glue_points)
        print("pushing points through generic jacobian isogenies:", t_jac_points)
        """
	
        return isogeny_chain


def evaluate_isogeny_jacobian(P,isogeny_chain):
    gluing = isogeny_chain[0]
    P = apply_gluing(P,gluing)
    for phi in isogeny_chain[1:-1]:
        P = apply_jacobian_isogeny(P,phi)
    splitting = isogeny_chain[-1]
    P = apply_splitting(P, splitting)
    return P
