
#cd C:\Users\BintuLabUser\Scope3AnalysisScripts\MERFISH_spot_analysis\smFISHAnalysisHuman && activate cellpose && python worker_smFISH__DCBB_human__03_23_2023_Scope3Analysis.py "_C" 50 "H4_"
from multiprocessing import Pool, TimeoutError
import time,sys
import os,sys
from tqdm import tqdm
master_analysis_folder = r'C:\Users\BintuLabUser\Scope3AnalysisScripts\MERFISH_spot_analysis\smFISHAnalysisHuman' #Scope3Microscope
sys.path.append(master_analysis_folder)
from ioMicro import *
def f(set_ifov_iQs=None):
    if True:
    
        
        try:
            asm = analysis_smFISH(data_folders = [r'\\192.168.0.10\bbfishdc13\DCBB_human__03_23_2023\H*',
            r'\\192.168.0.10\bbfishdc13\DCBB_human__03_23_2023\P*'],
                     save_folder =r'\\192.168.0.10\bbfishdc13\DCBB_human__03_23_2023_Analysis',
                     H0folder=  r'\\192.168.0.10\bbfishdc13\DCBB_human__03_23_2023\H0_*',exclude_H0=True)
            if set_ifov_iQs is None:
                return asm
            set_,ifov,iQs = set_ifov_iQs
            asm.set_set(set_)
            asm.set_fov(ifov)
            if False:
                final_segmentation(asm.fl_bk,
                                analysis_folder=asm.save_folder,
                                plt_val=True,
                                rescz = 3,trimz=2, resc=4,p99=8000)
            #iQ = 12
            #dic_psf = pickle.load(open(r'Y:\Glass_MERFISH\psf_Scope1_60x_cy5_clean.pkl','rb'))
            Qfolders = [qfld for qfld in asm.Qfolders if asm.set_ in os.path.basename(qfld) if os.path.isdir(qfld)]
            print(asm.Qfolders)
            if iQs is None: 
                iQs = np.arange(len(Qfolders))
            elif type(iQs) is str:
                tag = iQs
                iQs = [iqfld for iqfld,qfld in enumerate(Qfolders) if asm.set_ in os.path.basename(qfld)
                    if tag in os.path.basename(qfld) if os.path.isdir(qfld)]
            
            for iQ in iQs:
                asm.set_hybe(iQ)
                completed = asm.check_finished_file()
                if not completed:
                    
                    asm.get_background(force=False)
                    try:
                        asm.get_signal()
                        asm.compute_drift(sz=300)
                        asm.get_aligned_ims()
                        #imf = asm.im_sig__/asm.im_bk__
                        resc=5
                        asm.get_Xh(th = 750,s=30,delta_fit=4,dic_psf=None,subtract_bk=True,trim0=True,fr=None)
                        #get_Xh(self,th = 3,s=30,dic_psf=None,normalized=False,subtract_bk=False,trim0=True)
                        asm.dic_th = {0:4000,1:4000,2:4000}
                        asm.save_fits(icols=None,plt_val=True,save_max=True)
                    except:
                        print("Failed hybe",iQ)
        except:
            print("Failed")
    
    return None
if __name__ == '__main__':
    n = len(sys.argv)
    if n>=4:
        set_ = str(sys.argv[1])
        ifov = int(sys.argv[2])
        iQs = str(sys.argv[3])
        if iQs[0]=='[' and iQs[-1]==']': iQs=eval(iQs)
        set_ifov_iQs = (set_,ifov,iQs)
        f(set_ifov_iQs)
    elif n==1:
        asm = f(set_ifov_iQs=None)
        sets_sel = ['_C','_D']
        sets_ = np.array(['_'+os.path.dirname(fl).split('_')[-1]for fl in asm.fls_bk])
        #sets_,ifovs = np.unique(sets_,return_counts=True)
        ifovs = [np.sum(sets_==set_) for set_ in sets_sel]
        sets_ = sets_sel
        
        items = [(set_,ifov,None)for set_,nfov in zip(sets_,ifovs)
                                for ifov in range(nfov)]
        print(items)
        #f(items[40])
        if True:
            with Pool(processes=7) as pool:
                print('starting pool')
                start = time.time()
                result = pool.map(f, items)
                end = time.time()
            print(result,end-start)