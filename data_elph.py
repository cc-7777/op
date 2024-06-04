save_path='/home/publics/WORKSPACE/hongxi/2024/test/simple/com_nscf/bn.save/SAVE'
dipoles_path='/home/publics/WORKSPACE/hongxi/2024/test/simple/com_nscf/bn.save/2D_WR_WC'
bse_path ='/home/publics/WORKSPACE/hongxi/2024/test/simple/com_nscf/bn.save/2D_WR_WC'
ph_path='/home/publics/WORKSPACE/hongxi/2024/test/simple/dvscf/bn.save/SAVE'
from yambopy import *
import pickle

def save_phononenergies():
        """
        Read phonon frequencies in eV
        """
        ph_energies  = np.zeros([2304,12])
        
        for iq in range(2304):
            fil = ph_path+'/ndb.elph_gkkp_expanded_fragment_' + "%d"%(iq+1)
            database = Dataset(fil)
            ph_energies[iq] = np.sqrt(database.variables['PH_FREQS%d'%(iq+1)][:])*ha2ev
            database.close()
        with open('data_phon_enengy_bulk.pickle', 'wb') as f:
            pickle.dump(ph_energies,f)
def save_gkkp():

    gkkp = np.zeros([2304,12,16,16],dtype=np.complex64)

    for iq in range(2304):
        fil = ph_path+'/ndb.elph_gkkp_expanded_fragment_'+ "%d"%(iq+1)
        database = Dataset(fil)
        gkkp = database.variables['%s%d'%('ELPH_GKKP_Q',iq+1)][:]
        gkkp = np.swapaxes(gkkp[:,:,:,:,0] + I*gkkp[:,:,:,:,1],-1,1)
        database.close()
        with open('data_yelph_bulk_fragment_%d.pickle'%(iq+1), 'wb') as f:
            pickle.dump(gkkp,f)


if __name__ == "__main__":
    
    ylat = YamboLatticeDB.from_db_file(filename=save_path+'/ns.db1')
    
    #所需的式expanded_gkkp_01
    
    """
    with open('data_ylat_bulk.pickle', 'wb') as f:
        pickle.dump(ylat, f)'''
    
    '''ylat_noexp = YamboLatticeDB.from_db_file(filename=save_path+'/SAVE/ns.db1',Expand=False)
    with open('data_ylat_noexp_bulk.pickle', 'wb') as f:
        pickle.dump(ylat_noexp, f)
    yelph = YamboElectronPhononDB(ylat,folder_gkkp=save_path+'/SAVE',save=save_path+'/SAVE')
    """
    #with open('data_yelph_bulk.pickle', 'wb') as f:
        #pickle.dump(yelph, f)
    save_gkkp()
    """
    ydip = YamboDipolesDB(ylat,save=dipoles_path,filename='ndb.dipoles')
    with open('data_ydip_bulk.pickle', 'wb') as f:
        pickle.dump(ydip, f)
        save_phononenergies()
    """
    """    
    iQ=1
    yexc = YamboExcitonDB.from_db_file(ylat,filename='G:\Documents'+'/ndb.BS_diago_Q%d'%iQ)
    
    def get_excitons():
        #in eV energies[iQ][i_m]
        energies = np.zeros((1,9216))
        for iQ in range(1):
            yexc = YamboExcitonDB.from_db_file(ylat,filename=bse_path+'/ndb.BS_diago_Q%d'%(iQ+1)) 
           
            energies[iQ]    = yexc.eigenvalues.real 
        
        return energies
    def getA_T():
        A_exc=[]
        T_exc=[]
        for i in range(1,2):
            f=Dataset('G:\Documents'+'/ndb.BS_diago_Q'+str(i),'r') 
            BS_EIGENSTATES=f.variables['BS_EIGENSTATES'] # eigen-states Exciton eigenvectors are arranged as eigenvectors[i_exc, i_kvc]
            A_exc.append(np.array(BS_EIGENSTATES)[:,0:9216,0] + 1j*np.array(BS_EIGENSTATES)[:,0:9216,1])
            BS_TABLE=f.variables['BS_TABLE'] # BS_TABLE
            T_exc.append(np.array(BS_TABLE))  #T first index present Q,and t index k v c s_v s_c，last two is about spin.
    
        return A_exc,T_exc
    energy = get_excitons()
    A ,T = getA_T()
    with open('data_yexc_enengy_bulk.pickle_%d'%iQ, 'wb') as f:
            pickle.dump(energy,f)
    with open('data_yexc_A_bulk.pickle_%d'%iQ, 'wb') as f:
            pickle.dump(A,f)
    with open('data_yexc_T_bulk.pickle_%d'%iQ, 'wb') as f:
            pickle.dump(T,f)
    """
        
