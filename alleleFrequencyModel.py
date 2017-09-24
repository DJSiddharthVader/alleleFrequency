import numpy as np
import matplotlib.pyplot as plt

def af_sim(p,q,gens,pfit,qfit,pqfit):
    '''runs a simulation of allele frequency changes in a population over a number of generations, using initial allele
    frequencies and the relative fitness of a genotype. It assumes that population is constant, that the relative fitness
    of each genotype is constant and simple dominant/recessive gene interactions.
    :param p: Dominant allele frequency
    :param q: recessive allele frequency
    :param gens: number of generations
    :param pfit: dominant homozygous relative fitness
    :param qfit: recessive homozygous relative fitness
    :param pqfit:heterozygous relative fitness
    :return:tuple of allele frequencies at each generation
    '''
    p, q = round(p, 5), round(q, 5)
    if p+q != 1.0:
        raise ValueError
    pfreqs, qfreqs = [], []
    pfreqs.append(p)
    qfreqs.append(q)

    for x in range(gens):
        psquared, qsquared,pq2 = round(p*p,7), round(q*q,7), round(2 * p * q, 7)

        if round(psquared+qsquared+pq2,2) != 1.0:
            raise ValueError
        psquarednew = pfit*psquared
        pq2new = pqfit*pq2
        qsquarednew = qfit*qsquared
        totalafterselection = psquarednew+pq2new+qsquarednew
        psquarednew, pq2new, qsquarednew = psquarednew/totalafterselection, pq2new/totalafterselection, qsquarednew/totalafterselection
        p, q = psquarednew + pq2new/2, qsquarednew + pq2new/2
        pfreqs.append(p)
        qfreqs.append(q)
    return pfreqs, qfreqs

def AFSIM(p =0.01, q = 0.99, gens = 2000,pfit = 1.0, qfit = 0.996, pqfit = 0.998):
    '''using af_sim() it runs the allele frequency simulation varying the fitness levels for each genotype from 0.0-1.0
    with a step size of 0.1, while keeping the other two fitness levels constant.
    :param p: Dominant allele frequency
    :param q: recessive allele frequency
    :param gens: number of generations
    :param pfit: dominant homozygous relative fitness
    :param qfit: recessive homozygous relative fitness
    :param pqfit:heterozygous relative fitness
    :return:three graphs of the allele frequencies over time at varied fitness levels for each genotype
    '''
    gensarr = np.arange(gens)
    pfitl, qfitl , pqfitl = 0, 0, 0
    fits = [pfitl, pqfitl ,qfitl]
    W,w,E,e, R, r = [], [], [], [], [], []
    lists = [W,E,R,w,e,r]
    for d in range(3):
        a = fits[d]
        V, v = lists[d], lists[d+3]
        if d == 0:
            while a < 1.0:
                s, d = af_sim(p,q,gens,a,qfit,pqfit)
                V.append(s)
                v.append(d)
                a += 0.1
        elif d == 1:
            while a < 1.0:
                s, d = af_sim(p,q,gens,pfit,qfit,a)
                V.append(s)
                v.append(d)
                a += 0.1
        else:
            while a < 1.0:
                s, d = af_sim(p,q,gens,pfit,a,pqfit)
                V.append(s)
                v.append(d)
                a += 0.1
    genslist = np.arange(0,gens+1)

    f, (ax1,ax2) = plt.subplots(2)
    ax1.plot(genslist,lists[0][0], genslist, lists[0][1],genslist,lists[0][2], genslist,lists[0][3], genslist,lists[0][4], genslist,lists[0][5], genslist,lists[0][6], genslist,lists[0][7], genslist,lists[0][8], genslist,lists[0][9])
    ax2.plot(genslist,lists[3][0], genslist, lists[3][1],genslist,lists[3][2], genslist,lists[3][3], genslist,lists[3][4], genslist,lists[3][5], genslist,lists[3][6], genslist,lists[3][7], genslist,lists[3][8], genslist,lists[3][9])
    ax1.set_title("Dominant Allele Frequency Over Time Varying Dominant Homozygous Fitness", y=1.08)
    ax2.set_title("Recessive Allele Frequency Over Time Varying Dominant Homozygous Fitness", y=1.08)
    ax1.set_xlabel("Time (Number of Generations)")
    ax2.set_xlabel("Time (Number of Generations)")
    ax1.set_ylabel("Allele Frequency")
    ax2.set_ylabel("Allele Frequency")
    ax1.set_yticks([0.004, 0.008, 0.012, 0.016, 0.020])

    f1, (ax3, ax4) = plt.subplots(2)
    ax3.plot(genslist,lists[1][0], genslist, lists[1][1],genslist,lists[1][2], genslist,lists[1][3], genslist,lists[1][4], genslist,lists[1][5], genslist,lists[1][6], genslist,lists[1][7], genslist,lists[1][8], genslist,lists[1][9])
    ax4.plot(genslist,lists[4][0], genslist, lists[4][1],genslist,lists[4][2], genslist,lists[4][3], genslist,lists[4][4], genslist,lists[4][5], genslist,lists[4][6], genslist,lists[4][7], genslist,lists[4][8], genslist,lists[4][9])
    ax3.set_title("Dominant Allele Frequency Over Time Varying Heterozygous Fitness")
    ax4.set_title("Recessive Allele Frequency Over Time Varying Heterozygous Fitness")


    # f2, (ax5, ax6) = plt.subplots(2)
    # ax5.plot(genslist,lists[2][0], genslist, lists[2][1],genslist,lists[2][2], genslist,lists[2][3], genslist,lists[2][4], genslist,lists[2][5], genslist,lists[2][6], genslist,lists[2][7], genslist,lists[2][8], genslist,lists[2][9])
    # ax6.plot(genslist,lists[5][0], genslist, lists[5][1],genslist,lists[5][2], genslist,lists[5][3], genslist,lists[5][4], genslist,lists[5][5], genslist,lists[5][6], genslist,lists[5][7], genslist,lists[5][8], genslist,lists[5][9])
    # ax5.set_title("Dominant Allele Frequency Over Time Varying Recessive Homozygous Fitness")
    # ax6.set_title("Recessive Allele Frequency Over Time Varying Recessive Homozygous Fitness")
    #
    # plt.legend(("Relative Varied Genotype Fitness 0.1","Relative Varied Genotype Fitness 0.2","Relative Varied Genotype Fitness 0.3","Relative Varied Genotype Fitness 0.4","Relative Varied Genotype Fitness 0.5","Relative Varied Genotype Fitness 0.6","Relative Varied Genotype Fitness 0.7","Relative Varied Genotype Fitness 0.8","Relative Varied Genotype Fitness 0.9","Relative Varied Genotype Fitness 1.0"))

    plt.legend(("Homozygous Dominant Fitness 0.1","Homozygous Dominant Fitness 0.2","Homozygous Dominant Fitness  0.3","Homozygous Dominant Fitness  0.4","Homozygous Dominant Fitness  0.5","Homozygous Dominant Fitness  0.6","Homozygous Dominant Fitness  0.7","Homozygous Dominant Fitness  0.8","Homozygous Dominant Fitness  0.9","Homozygous Dominant Fitness 1.0"))
    plt.show()

AFSIM()

def af_simgraph(p,q,gens,pfit,pqfit, qfit, title):
    '''runs a simulation of allele frequency changes in a population over a number of generations, using initial allele
    frequencies and the relative fitness of a genotype. It assumes that the relative fitness
    of each genotype is constant and simple dominant/recessive gene interactions.
    :param p: Dominant allele frequency
    :param q: recessive allele frequency
    :param gens: number of generations
    :param pfit: dominant homozygous relative fitness
    :param qfit: recessive homozygous relative fitness
    :param pqfit:heterozygous relative fitness
    :return:tuple of allele frequencies at each generation
    '''
    p, q = round(p, 5), round(q, 5)
    if p+q != 1.0:
        raise ValueError
    pfreqs, qfreqs = [], []
    pfreqs.append(p)
    qfreqs.append(q)

    for x in range(gens):
        psquared, qsquared,pq2 = round(p*p,7), round(q*q,7), round(2 * p * q, 7)

        if round(psquared+qsquared+pq2,2) != 1.0:
            raise ValueError
        psquarednew = pfit*psquared
        pq2new = pqfit*pq2
        qsquarednew = qfit*qsquared
        totalafterselection = psquarednew+pq2new+qsquarednew
        psquarednew, pq2new, qsquarednew = psquarednew/totalafterselection, pq2new/totalafterselection, qsquarednew/totalafterselection
        p, q = psquarednew + pq2new/2, qsquarednew + pq2new/2
        pfreqs.append(p)
        qfreqs.append(q)
    genslist = np.arange(0,gens+1)
    plt.plot(genslist, pfreqs)
    plt.plot(genslist, qfreqs)
    plt.suptitle(title)
    plt.xlabel("Time (# of generations)")
    plt.ylabel("Allele Frequency")
    plt.legend(("Dominant Allele", "Recessive Allele"))
    plt.axis([0,gens,-0.1,1.1])
    plt.show()

#af_simgraph(0.51,0.49,1000,1,0.99,1.0, "Allele Frequencies for Underdominance")