from multiprocessing import Pool, TimeoutError
import time,sys
import os,sys,numpy as np
#sys.path.append(r'C:\Users\BintuLab\Dropbox\MERFISH_DC_SCOPE3')
master_analysis_folder = r'C:\Users\BintuLab\Scope4AnalysisScripts\MERFISH_Spot_Analysis\Analysis_1500gns_MERFISH'
sys.path.append(master_analysis_folder)
from ioMicro import *
def compute_drift(save_folder,fov,all_flds,set_,redo=False):
    drift_fl = save_folder+os.sep+'drift_'+fov.split('.')[0]+'--'+set_+'.pkl'
    if not os.path.exists(drift_fl) or redo:
        ims = [read_im(fld+os.sep+fov) for fld in all_flds]
        ncols,sz,sx,sy = ims[0].shape
        im_ref = np.array(ims[len(ims)//2][-1],dtype=np.float32)
        drifts = [get_txyz(im[-1],im_ref,sz_norm=20, sz=500) for im in tqdm(ims)]
        
        pickle.dump([drifts,all_flds,fov],open(drift_fl,'wb'))
def compute_fits(save_folder,fov,all_flds,redo=False):
    for fld in tqdm(all_flds):
        im_ = read_im(fld+os.sep+fov)
        ncols = len(im_)
        for icol in range(ncols-1):
            #icol=2
            tag = os.path.basename(fld)
            save_fl = save_folder+os.sep+fov.split('.')[0]+'--'+tag+'--col'+str(icol)+'__Xhfits.npz'
            if not os.path.exists(save_fl) or redo:
                im__ = np.array(im_[icol],dtype=np.float32)
                im_n = norm_slice(im__,s=30)
                Xh = get_local_max(im_n,500,im_raw=im__,dic_psf=None,delta=1,delta_fit=3,dbscan=True,
                      return_centers=False,mins=None,sigmaZ=1,sigmaXY=1.5)
                np.savez_compressed(save_fl,Xh=Xh)
def main_f(set_ifov):
    save_folder =r'Y:\DCBBL1_3_2_2023\MERFISH_Analysis'
    if not os.path.exists(save_folder): os.makedirs(save_folder)
    all_flds = [r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\Controls\H0_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\Controls\H1_Igfbp_Aldh1l1_Ptbp1_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H1_MER_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H2_MER_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H3_MER_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H4_MER_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H5_MER_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H6_MER_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H7_MER_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H8_MER_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H9_MER_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H10_MER_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H11_MER_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H12_MER_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H13_MER_set1',
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H14_MER_set1'
       r'\\192.168.0.3\bbfishdc10_2\DCBBL1_3_2_2023\H15_MER_set1'
       r'\\192.168.0.3\bbfishdc10_2\DCBBL1_3_2_2023\H16_MER_set1']
    set_,ifov = set_ifov
    all_flds = [fld.replace('_set1',set_) for fld in all_flds]
    fovs_fl = save_folder+os.sep+'fovs__'+set_+'.npy'
    if not os.path.exists(fovs_fl):
        fls = glob.glob(all_flds[0]+os.sep+'*.zarr')
        fovs = [os.path.basename(fl) for fl in fls]
        np.save(fovs_fl,fovs)
    else:
        fovs = np.load(fovs_fl)
    if ifov<len(fovs):
        fov = fovs[ifov]
        print("Computing fitting on: "+str(fov))
        compute_fits(save_folder,fov,all_flds,redo=False)
        print("Computing drift on: "+str(fov))
        compute_drift(save_folder,fov,all_flds,set_,redo=True)
    return set_ifov
if __name__ == '__main__':
    # start 4 worker processes
    items = [(set_,ifov)for set_ in ['_set1','_set2','_set3','_set4']
                        for ifov in range(1000)]
                        
    #main_f(['_set1',50000])
    if True:
        with Pool(processes=35) as pool:
            print('starting pool')
            result = pool.map(main_f, items)
            if False:
                start = time.time()
                result = pool.map(main_f, items)
                end = time.time()
                print('Sleeping for a bit')
                print(result,end-start)
                time.sleep(600)
    