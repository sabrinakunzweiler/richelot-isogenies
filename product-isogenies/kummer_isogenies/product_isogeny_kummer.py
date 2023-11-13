import time

from sage.all import ZZ, inverse_mod
#from utilities.strategy import optimised_strategy
from kummer_isogenies.formulae_kummer import *
from kummer_isogenies.doubling_kummer import *
from kummer_isogenies.gluing_isogeny_kummer import *
from kummer_isogenies.generic_isogeny_kummer import *
from kummer_isogenies.splitting_isogeny_kummer import *



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

def isogeny_chain_kummer(kernel, n):
        """
        We compute a (2^n,2^n)-chain.
        The input is expected to be a (2^(n+1),2^(n+1))-kernel.
        """

        # Extract points from kernel
        (P1, Q1, P2, Q2) = kernel

        Tp1 = [P1, P2]
        Tp2 = [Q1, Q2]

        # Store chain of (2,2)-isogenies
        # we store matrices and transformations
        isogeny_chain = []
        phi = 0

        # Bookkeeping for optimal strategy
        strat_idx = 0
        level = [0]
        ker = [Tp1, Tp2]
        kernel_elements = [ker]
        strategy = optimised_strategy_kummer(n, mul_c=1.4)
        t_double_ell = 0
        t_double_kum = 0
        t_glue_points = 0
        t_kum_points = 0
        t_compute_glue = 0
        t_compute_kummer_constants = 0
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
                phi = gluing_kummer_isogeny(ker)
                t_compute_glue = t_compute_glue + time.process_time() -t0
            elif k == 1:
                #in this case the usual doubling strategy doesn't work
                t0 = time.process_time()
                phi = first_kummer_isogeny(ker, isogeny_chain[-1], ker2)
                t_compute_kummer_constants = t_compute_kummer_constants + time.process_time() -t0
            elif k == n-1:
                #is this always 4-torsion?
                t0 = time.process_time()
                phi = splitting_kummer_isogeny(ker, phi)
                #isogeny_chain.append(splitting_isogeny)
                t_compute_split = t_compute_split + time.process_time() -t0
            else:
                t_double_kum = t_double_kum +  time.process_time() - t0
                t0 = time.process_time()
                phi = kummer_isogeny(ker, isogeny_chain[-1])
                t_compute_kummer_constants = t_compute_kummer_constants + time.process_time() -t0
                        
            # Update the chain of isogenies
            isogeny_chain.append(phi)

            # Remove elements from list
            kernel_elements.pop()
            level.pop()

            # Push through points for the next step
            t0 = time.process_time()
            if k == 0:
                kernel_elements = [[apply_gluing(ker[0],phi), apply_gluing(ker[1],phi)] for ker in kernel_elements]
                #for the next isogeny we need to precompute kernel2
                ker2 = [apply_gluing(Tp1,phi), apply_gluing(Tp2,phi)]
                t_glue_points = t_glue_points + time.process_time() - t0
            elif k != n-1:
                kernel_elements = [(apply_kummer_isogeny(ker[0],phi), apply_kummer_isogeny(ker[1],phi)) for ker in kernel_elements]
                t_kum_points = t_kum_points + time.process_time() - t0
        """
        print("doublings on elliptic curve:", t_double_ell)
        print("doublings on kummer surface:", t_double_kum)
        print("compute the gluing constants:", t_compute_glue)
        print("compute the kummer constants:", t_compute_kummer_constants)
        print("compute the splitting constants:", t_compute_split)
        print("pushing points through gluing:", t_glue_points)
        print("pushing points through generic kummer isogenies:", t_kum_points)
        """
	
        return isogeny_chain


def evaluate_isogeny_kummer(P,isogeny_chain):
        gluing = isogeny_chain[0]
        P = apply_gluing(P,gluing)
        for phi in isogeny_chain[1:-1]:
            P = apply_kummer_isogeny(P,phi)
        splitting = isogeny_chain[-1]
        P = apply_splitting(P, splitting)
        return P
