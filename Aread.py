from yambopy import *
import numpy as np
import matplotlib.pyplot as plt
import h5py
import pickle
save_path='/home/publics/WORKSPACE/hongxi/2024/test/simple/com_nscf/bn.save/SAVE'
bse_path ='/home/publics/WORKSPACE/hongxi/2024/test/simple/com_nscf/bn.save/2D_WR_WC'
#elph_path='/home/cc/databases_yambopy/ELPH_saves/YAMBO_saves'
#dipoles_path='/home/cc/databases_yambopy/BSE_saves/BSE_databases'
    
def foldstd(klist):  
    #when x,y,z out of box ,get return
    ylpoxs=np.zeros(2304)
    ylpoys=np.zeros(2304)
    ylpozs=np.zeros(2304)
    for i in range(2304):
        ylpoxs[i]=klist[i][0]
        ylpoys[i]=klist[i][1]
        ylpozs[i]=klist[i][2]
        while ylpoxs[i] > 23 or ylpoxs[i] < 0:
            if ylpoxs[i] >23:
                ylpoxs[i]=ylpoxs[i]-24
            else:
                ylpoxs[i]=ylpoxs[i]+24
        while ylpoys[i] > 23 or ylpoys[i] < 0:
            if ylpoys[i] >23:
                ylpoys[i]=ylpoys[i]-24
            else:
                ylpoys[i]=ylpoys[i]+24
        while ylpozs[i] > 3 or ylpozs[i] < 0:
            if ylpozs[i] >3:
                ylpozs[i]=ylpozs[i]-4
            else:
                ylpozs[i]=ylpozs[i]+4
    newlist=np.zeros((2304,3))
    for i in range(2304):
        newlist[i]=[ylpoxs[i],ylpoys[i],ylpozs[i]]
    return newlist

def kmatch(k,map):
    '''change the three dimension to two dimension,fold the k into a standard square.k use the 
    reduced,and multiply by 24.'''
    kx = round(k[0]*24)
    ky = round(k[1]*24)
    kz = round(k[2]*4)
    while kx > 23 or kx < 0:
            if kx >23:
                kx=kx-24
            else:
                kx=kx+24
    while ky > 23 or ky < 0:
            if ky >23:
                ky=ky-24
            else:
                ky=ky+24
    while kz > 3 or kz < 0:
            if kz >3:
                kz=kz-4
            else:
                kz=kz+4
    ikorigin=map[kz*24*24+ky*24+kx]  #map is the b to origin.
    return ikorigin

def getchanged_ylk():  #get the ylat kpoints red_kpoints multiply by 12 and round.
    ylkpoints=np.zeros((2304,3))
    for i in range(2304):
        ylkpoints[i][0]=round(ylat.red_kpoints[i][0]*24)
        ylkpoints[i][1]=round(ylat.red_kpoints[i][1]*24)
        ylkpoints[i][2]=round(ylat.red_kpoints[i][2]*4)
    '''ylpoxs=np.zeros((144)) 
    ylpoys=np.zeros((144))
    for i in range(144):
        ylpoxs[i]=ylkpoints[i][0]
    for i in range(144):
        ylpoys[i]=ylkpoints[i][1]
    plt.figure(figsize=(10,10),dpi=100)
    plt.scatter(ylpoxs,ylpoys)
    plt.show()'''
    #print(ylpoxs)
    return ylkpoints

def getexciqpoints():
    #exciqponits=[] 
    #直接从r_setup读    
    '''for i in range(1,184):#print(yexc)
        yexc = YamboExcitonDB.from_db_file(ylat,filename=bse_path+'/ndb.BS_diago_Q%d'%i)
        if yexc.car_qpoint is  None:
            exciqponits.append([0,0,0])
        else :
            exciqponits.append(yexc.car_qpoint*yexc.lattice.alat)'''
    file = open('/home/publics/WORKSPACE/hongxi/2024/test/simple/data_qlist.pickle','rb')  # 以二进制读模式（rb）打开pkl文件
    exciqponits = pickle.load(file)
    return exciqponits

def getstandardklist():
     #creat a standard square 12x12 from zero to 11
    klist=np.zeros((2304,3))
    for k in range(4):
        for j in range(24):
            for i in range(24):
                klist[k*24*24+24*j+i]=[i,j,k]
    #print(klist)
    return klist

def getA_T():
    A_exc=[]
    T_exc=[]
    for i in range(1,184):
        f=Dataset(bse_path+'/ndb.BS_diago_Q'+str(i),'r') 
        BS_EIGENSTATES=f.variables['BS_EIGENSTATES'] # eigen-states Exciton eigenvectors are arranged as eigenvectors[i_exc, i_kvc]
        A_exc.append(np.array(BS_EIGENSTATES)[0:10,0:9216,0] + 1j*np.array(BS_EIGENSTATES)[0:10,0:9216,1])
        BS_TABLE=f.variables['BS_TABLE'] # BS_TABLE
        T_exc.append(np.array(BS_TABLE))  #T first index present Q,and t index k v c s_v s_c，last two is about spin.
    filename='A_exc'
    dataname='A_exc_data'
    with h5py.File("%s.hdf5"%filename, 'w') as f:  # 写入的时候是‘w’
        f.create_dataset("%s"%dataname, data=A_exc, compression="gzip", compression_opts=5)
    return A_exc,T_exc

def point_matching2(a,b,double_check=True,debug=False,eps=0.001):
    """
    Matches the points of list a to the points of list b
    using a nearest neighbour finding algorithm

    Arguments:

        double_check: after the nearest neighbours are assigned check further
        if the distance between points is within the precision eps

        eps: precision for the double check (default: 1e-8)

    """
    #karma
    from scipy.spatial import cKDTree
    from time import time
    a = np.array(a)
    b = np.array(b)
    start_time = time()

    #initialize thd kdtree
    kdtree = cKDTree(a, leafsize=10)
    map_b_to_a = []
    for xb in b:
        current_dist,index = kdtree.query(xb, k=1, distance_upper_bound=6)
        map_b_to_a.append(index)
    map_b_to_a = np.array(map_b_to_a)
    
    if debug:
        print("took %4.2lfs"%(time()-start_time))

    if double_check:
        for ib,ia in enumerate(map_b_to_a):
            dist = np.linalg.norm(a[ia]-b[ib])
            if dist > eps:
                raise ValueError('point a %d: %s is far away from points b %d: %s  dist: %lf'%(ia,str(a[ia]),ib,str(b[ib]),dist))

    return map_b_to_a

def ikvc_to_it(ik,vk,ck,table):
        for it, t in enumerate(table):
            if t[0]-1==ik and t[1]-7==vk and t[2]-7==ck: return it
        

def it_to_ikvc(it,table):#输出是个元组,分别是ik,iv,ic，也就是k点序数，价带序数和导带序数
    t = table[it]
    return t[0]-1,t[1]-7,t[2]-7


def Q_expand_to_noexp(iQ,list1):
        
    return(list1[iQ])

def get_A_multi1(im,i_n,iQ,iq,ik,iv,ic,ic2):
    it = ikvc_to_it(ik,iv,ic,table)
    it2 = ikvc_to_it(ik,iv,ic2,table)
    Qk=exciqponits[iQ]
    qk=ylat.red_kpoints[iq]
    iQplusq=kmatch(Qk+qk,map_fold_tostd)
    noexp_iQplusq = Q_expand_to_noexp(iQplusq,list1) 
    A_dia_Sm_Qplusq_vk_c_kplusqplusQ =  np.conj(A_exc[noexp_iQplusq][im][it])
    A_Sn_Q_vk_c2_kplusQ = A_exc[iQ][i_n][it2]
    A_result = A_dia_Sm_Qplusq_vk_c_kplusqplusQ*A_Sn_Q_vk_c2_kplusQ
    return A_result

def get_A_multi2(im,i_n,iQ,iq,ik,iv,iv2,ic):
    it2 = ikvc_to_it(ik,iv2,ic,table)
    Qk=exciqponits[iQ]
    qk=ylat.red_kpoints[iq]
    k_k = ylat.red_kpoints[ik]
    iQplusq=kmatch(Qk+qk,map_fold_tostd)
    noexp_iQplusq = Q_expand_to_noexp(iQplusq,list1) 
    i_kminq = kmatch(k_k-qk,map_fold_tostd)
    it = ikvc_to_it(i_kminq,iv,ic,table)
    A_dia_Sm_Qplusq_vk_c_kplusqplusQ =np.conj(A_exc[noexp_iQplusq][im][it])
    A_Sn_Q_v2k_c_kplusQ = A_exc[iQ][i_n][it2]
    A_result = A_dia_Sm_Qplusq_vk_c_kplusqplusQ*A_Sn_Q_v2k_c_kplusQ
    return A_result

def get_G(i_m,i_n,iQ,iq,i_nu):
    multi1 = 0
    multi2 = 0
    multi3 = 0
    for ik in range(2304):
        Qk=exciqponits[iQ]
        k_k=ylat.red_kpoints[ik]
        iQplusk=kmatch(Qk+k_k,map_fold_tostd)
        qk=ylat.red_kpoints[iq]
        i_kminq = kmatch(k_k-qk,map_fold_tostd)

        for iv in range(0,2):
            for ic in range(2,4):
                for ic2 in range(2,4):
                    A1 = get_A_multi1(i_m,i_n,iQ,iq,ik,iv,ic,ic2)
                    g1 = gkkp[iq][iQplusk][i_nu][ic][ic2]
                    multi1 = A1*g1 + multi1 
        for iv in range(0,2):
            for iv2 in range(0,2):
                for ic in range(2,4):
                    A2 = get_A_multi2(i_m,i_n,iQ,iq,ik,iv,iv2,ic2)
                    g2 = gkkp[iq][i_kminq][i_nu][iv][iv2]
                    multi2 = A2*g2 + multi2
        multi3 = multi3+multi1-multi2
        
    print(multi3)
    return multi3



def get_gkkp():  #缩减到只有4个带，也就是2304*2304*12*4*4
    gkkp = np.zeros((2304,2304,12,4,4),dtype=np.complex64)
    for iq in range(2304):
        file = open('/home/publics/WORKSPACE/hongxi/2024/test/simple/elph_data/data_yelph_bulk_fragment_%d.pickle'%(iq+1),'rb')  # 以二进制读模式（rb）打开pkl文件
        data = pickle.load(file)
        for i_m in range(4):
            for i_n in range(4):
                for i_k in range(2304):
                    for i_nu in range(12):
                        gkkp[iq][i_k][i_nu][i_m][i_n] = data[i_k][i_nu][i_m+6][i_n+6]
    filename = 'gkkp'
    dataname = 'gkkp_data'

    with h5py.File("%s.hdf5"%filename, 'w') as f:  # 写入的时候是‘w’
            f.create_dataset("%s"%dataname, data=gkkp, compression="gzip", compression_opts=5)
  
                
        
    return gkkp

def read_gkkp(filename,dataname):
    with h5py.File("%s.hdf5"%filename, 'r') as f:  # 读取的时候是‘r’
        gkkp = f.get("%s"%dataname)[:]
    
    return gkkp

def read_A(filename,dataname):
    with h5py.File("%s.hdf5"%filename, 'r') as f:  # 读取的时候是‘r’
        A_T = f.get("%s"%dataname)[:]
    
    return A_T


def get_G1(i_m,i_n,i_Q,i_q,i_nu):
    multi1 = 0
    multi2 = 0
    multi3 = 0
    for ik in range(1):
        Qk=exciqponits[i_Q]
        k_k=ylat.red_kpoints[ik]
        iQplusk=kmatch(Qk+k_k,map_fold_tostd)
        qk=ylat.red_kpoints[i_q]
        i_kminq = kmatch(k_k-qk,map_fold_tostd)

        for iv in range(0,2):
            for ic in range(2,4):
                for ic2 in range(2,4):
                    A1 = get_A_multi1(i_m,i_n,i_Q,i_q,ik,iv,ic,ic2)
                    g1 = gkkp[i_q][iQplusk][i_nu][ic][ic2]
                    multi1 = A1*g1 + multi1
        for iv in range(0,2):
            for iv2 in range(0,2):
                for ic in range(2,4):
                    A2 = get_A_multi2(i_m,i_n,i_Q,i_q,ik,iv,iv2,ic2)
                    g2 = gkkp[i_q][i_kminq][i_nu][iv][iv2]
                    multi2 = A2*g2 + multi2
        multi3 = multi3+multi1-multi2
        
    print(multi3)
    return multi3









def get_excitons():
    #in eV energies[iQ][i_m]
    energies = np.zeros((183,9216))
    for iQ in range(183):
        yexc = YamboExcitonDB.from_db_file(ylat,filename=bse_path+'/ndb.BS_diago_Q%d'%(iQ+1)) 
           
        energies[iQ]    = yexc.eigenvalues.real 
        
    return energies


def save_phononenergies():
        """
        Read phonon frequencies in eV
        """
        ph_energies  = np.zeros([2304,12])
        
        for iq in range(2304):
            fil = '/home/publics/WORKSPACE/hongxi/2024/test/simple/dvscf/bn.save/SAVE'+'/ndb.elph_gkkp_expanded_fragment_' + "%d"%(iq+1)
            database = Dataset(fil)
            ph_energies[iq] = np.sqrt(database.variables['PH_FREQS%d'%(iq+1)][:])*ha2ev
            database.close()
        with open('data_phon_enengy_bulk.pickle', 'wb') as f:
            pickle.dump(ph_energies,f)


def get_phonon():
    #yelph.ph_energies[iq][i_nu]   in eV
    file = open('/home/publics/WORKSPACE/hongxi/2024/test/simple/data_phon_enengy_bulk.pickle','rb')  # 以二进制读模式（rb）打开pkl文件
    ph_energies = pickle.load(file)

    return ph_energies

def Bosedis(energy,T,min_exp=-100.):
    k = 8.6173333e-05  
    argument=energy/(k*T)
    argument=np.max([argument,min_exp])
    N = 1/(np.exp(argument)-1)
    return N

def Boltz_exdis(e1,energy,T,min_exp=-100.):
    k = 8.6173333e-05
    argument=-(energy-e1)/(k*T)
    argument=np.max([argument,min_exp])
    N = np.exp(argument)
    return N


def delta(energy,broad = 4e-3):
    n = np.exp((-1)*energy**2/broad**2)
    m = n/(np.pi**(1/2)*broad)
    return m

def deltasp(energy,h=0.004):
    if energy >h:
        height = 1/h
    else:
        height = 0
    return height
    
def squarelimit(x,x0,s):
    if abs(x-x0)>s:
        h = x-x0
    else:
        h = s
    return 1/h**2
def abs2(x):
    return x.real**2 + x.imag**2

def gaussian(x,x0,s,max_exp=50.,min_exp=-100.):
    height=1./(np.sqrt(2.*np.pi)*s)
    argument=-0.5*((x-x0)/s)**2
    #Avoiding undeflow errors...
    argument=np.max([argument,min_exp])
    #when the number is too small,change it to min_exp.
    return height*np.exp(argument)

def read_G(filename,dataname):
    with h5py.File("%s.hdf5"%filename, 'r') as f:  # 读取的时候是‘r’
        G_data = f.get("%s"%dataname)[:]
    
    return G_data

def save_G(filename,dataname):
        """Compute G and dump it to file"""
        G_noexp = np.zeros([8,8,183,2304,12],dtype=np.complex64)
        for i_m in range(8):
            for i_n in range(8):
                for i_Q in range(183):
                    for i_q in range(2304):
                        for i_nu in range(12):
                            G_noexp[i_m,i_n,i_Q,i_q,i_nu]=get_G(i_m,i_n,i_Q,i_q,i_nu)
        
        

        
        with h5py.File("%s.hdf5"%filename, 'w') as f:  # 写入的时候是‘w’
            f.create_dataset("%s"%dataname, data=G_noexp, compression="gzip", compression_opts=5)
        return True


def lifetime(i_n,iQ,T):
    plus = 0
    for i_m in range(5):
        for i_nu in range(6):
            for iq in range(144):
                G = G_data[i_m][i_n][iQ][iq][i_nu]
                EnQ = excienergies[iQ][i_n]   #n,Q exciton energy
                e1 = excienergies[iQ][0]
                Qk=exciqponits[iQ]
                qk=ylat.red_kpoints[iq]
                iQplusq=kmatch(Qk+qk,map_fold_tostd)
                noexp_iQplusq = Q_expand_to_noexp(iQplusq,list1)
                EmQplusq = excienergies[noexp_iQplusq][i_m]
                Eph_en = phononenergies[iq][i_nu]
                N_nu_q = Bosedis(Eph_en,T)
                F_m_Qplusq = Boltz_exdis(e1,EmQplusq,T)
                plus += abs2(G)*((N_nu_q + 1 + F_m_Qplusq)*gaussian((EnQ-EmQplusq-Eph_en),0,9e-3)+(N_nu_q - F_m_Qplusq )*gaussian(EnQ-EmQplusq+Eph_en,0,4e-3))
    return 1/(((2*np.pi)/(hbar*144))*plus)

def light(energy,T):
    plus = 0
    for i_n in range(5):
        for i_nu in range(5,6):
            for i_m in range(5):
                for iQ in range(19):    
                    EnQ = excienergies[iQ][i_n]+3.72   #n,Q exciton energy
                    Qk=exciqponits[iQ]
                    qk=[Qk[0]*(-1),Qk[1]*(-1),Qk[2]*(-1)]
                    e1 = excienergies[iQ][0]+3.72
                    iq = kmatch(qk,map_fold_tostd)
                    G = G_data[i_m][i_n][iQ][iq][i_nu]
                    Em = excienergies[0][i_m]+3.72
                    Eph_en = phononenergies[iq][i_nu]
                    NnQ = Boltz_exdis(e1,EnQ,T)
                    Nph_n_u_Q = Bosedis(Eph_en,T)
                    i_k,i_v,i_c=it_to_ikvc(i_m,table)
                    dip_of_k = np.abs(ydip.dipoles[i_k,0,i_c,i_v])+np.abs(ydip.dipoles[i_k,1,i_c,i_v])
                    plus += np.square(dip_of_k)*abs2(G)*gaussian(energy+Eph_en-EnQ,0,9e-3)*NnQ*(1+Nph_n_u_Q)/squarelimit(energy,Em,0.004)
    return plus

def light2(energy,T):
    plus = 0
    for i_n in range(5):
        for i_nu in range(6):
            for i_m in range(5):
                for iQ in range(19):    
                    EnQ = excienergies[iQ][i_n]+3.72   #n,Q exciton energy
                    Qk=exciqponits[iQ]
                    qk=[Qk[0]*(-1),Qk[1]*(-1),Qk[2]*(-1)]
                    e1 = excienergies[iQ][0]+3.72
                    iq = kmatch(qk,map_fold_tostd)
                    G = G_data[i_m][i_n][iQ][iq][i_nu]
                    Em = excienergies[0][i_m]+3.72
                    Eph_en = phononenergies[iq][i_nu]
                    NnQ = Boltz_exdis(e1,EnQ,T)
                    Nph_n_u_Q = Bosedis(Eph_en,T)
                    Nh_w = Bosedis(energy,T)
                    Ne_m = Boltz_exdis(e1,Em,T)
                    i_k,i_v,i_c=it_to_ikvc(i_m)
                    dip_of_k = np.abs(ydip.dipoles[i_k,0,i_c,i_v])+np.abs(ydip.dipoles[i_k,1,i_c,i_v])
                    plus1 = Nh_w*gaussian(energy+Eph_en-EnQ,0,9e-3)/squarelimit(energy,Em,0.004)
                    plus2 = np.square(np.pi)*(Nph_n_u_Q-1)*gaussian(energy+Eph_en-EnQ,0,9e-3)*np.square(gaussian(EnQ-Em-Eph_en,0,4e-2))
                    plus3 = Ne_m*gaussian(EnQ-Eph_en-Em,0,9e-3)/squarelimit(energy,Em,0.004)
                    plus4 = Ne_m*(1+Ne_m)*gaussian(EnQ-Eph_en-Em,0,9e-3)/(Em-energy)
                    plus += np.square(dip_of_k)*abs2(G)*(plus1+plus2-plus3-plus4)
    return plus

                

if __name__ == "__main__":


    #Create "lattice" object by reading the ns.db1 database inside the yambo SAVE
    ylat = YamboLatticeDB.from_db_file(filename=save_path+'/ns.db1')
    #yelph = YamboElectronPhononDB(ylat,folder_gkkp=elph_path+'/SAVE',save=elph_path+'/SAVE')
    #gkkp改为读取pickle文件，gkkp[iq][iQplusk][i_nu][ic][ic2]，总共有2304*2304*12*16*16个数据
    #这里的导带为i_n=0到i_n=7，价带是i_n=8到i_n=15，这里我们只取两条导带和两条价带
    #gkkp = get_gkkp()
    gkkp = read_gkkp(filename='gkkp',dataname='gkkp_data') 
    exciqponits=getexciqpoints()     
    #A_exc,T_exc = getA_T()  #A[iQ][i_exc][i_kvc]
    A_exc = read_A(filename='A_exc',dataname='A_exc_data')
    ylkpoints=getchanged_ylk()
    standklist= getstandardklist()
    # fold ylatkpoints to standard,find the map to standard list.
    map_fold_tostd=point_matching2(foldstd(ylkpoints),standklist)
    #ydip = YamboDipolesDB(ylat,save=dipoles_path,filename='ndb.dipoles')
    list1 = ylat.kpoints_indexes 
    iQ1 = 0 #取名区分开
    yexc = YamboExcitonDB.from_db_file(ylat,filename=bse_path+'/ndb.BS_diago_Q%d'%(iQ1+1))
    table = yexc.table
        
    '''k1 = [0.16666,0,0]
    print(kmatch(k1,map_fold_tostd)) #show the k to the origin ik.k use the red_kpoints.'''
    # gkkp[q][k][mode][bnd1][bnd2],c (bnd1)to c'(bnd2)
    excienergies = get_excitons()
    #save_phononenergies()
    phononenergies = get_phonon()
    hbar = 6.582e-16 #h/2pi eV*s
    #get_G1(0,0,0,0,0)
    aa = save_G(filename='G_data',dataname='G_noexp_data')
    '''
    G_data = read_G(filename='G_data',dataname='G_noexp_data')
    a = [i/400+5.70 for i in range(200)]
    b = []
    for i in range(200):
        b.append(light2(a[i],200))
    
    output = open('data10.txt','w+')
    for j in range(200):
        output.write(str(a[j]))
        output.write(' ')
        output.write(str(b[j]))
        output.write('\n')
    output.close()
    plt.plot(a,b)
    plt.show
    #print(lifetime(0,4,300))'''
    
    

    
 
    
    
     
    
    
    
                       
