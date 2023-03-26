
#cd "C:\Users\BintuLab\Scope4AnalysisScripts\MERFISH_Spot_Analysis\Analysis_1500gns_MERFISH"&&activate cellpose&&python worker3.py

from multiprocessing import Pool, TimeoutError
import time,sys
import os,sys,numpy as np

#sys.path.append(r'C:\Users\BintuLab\Dropbox\MERFISH_DC_SCOPE3')
master_analysis_folder = r'C:\Users\BintuLab\Scope4AnalysisScripts\MERFISH_Spot_Analysis\Analysis_1500gns_MERFISH'
sys.path.append(master_analysis_folder)
from ioMicro import *


def compute_drift(save_folder,fov,all_flds,set_,redo=False):
    """
    save_folder where to save to - analysis_fodler
    fov - i.e. Conv_zscan_005.zarr
    all_flds - folders that contain eithger the MERFISH bits or control bits or sm
    
    """
    #print(len(all_flds))
    #print(all_flds)
    drift_fl = save_folder+os.sep+'drift_'+fov.split('.')[0]+'--'+set_+'.pkl'
    iiref = None
    previous_drift = {}
    if not os.path.exists(drift_fl):
        redo = True
    else:
        drifts_,all_flds_,fov_ = pickle.load(open(drift_fl,'rb'))
        all_tags_ = np.sort([os.path.dirname(fld)for fld in all_flds_])
        all_tags = np.sort([os.path.dirname(fld)for fld in all_flds])
        iiref = np.argmin([np.sum(np.abs(drift[0]))for drift in drifts_])
        previous_drift = {tag:drift for drift,tag in zip(drifts_,all_tags_)}
        if not (len(all_tags_)==len(all_tags)):
            redo = True
        else:
            if not np.all(all_tags_==all_tags):
                redo = True
    if redo:
        print("Computing drift...")
        ims = [read_im(fld+os.sep+fov) for fld in all_flds] #map the image
        ncols,sz,sx,sy = ims[0].shape
        if iiref is None: iiref = len(ims)//2
        im_ref = np.array(ims[iiref][-1],dtype=np.float32)
        all_tags = np.sort([os.path.dirname(fld)for fld in all_flds])
        drifts = [previous_drift.get(tag,get_txyz(im[-1],im_ref,sz_norm=20, sz=400)) 
                    for im,tag in zip(tqdm(ims),all_tags)]
        
        pickle.dump([drifts,all_flds,fov],open(drift_fl,'wb'))
def compute_fits(save_folder,fov,all_flds,redo=False,ncols=4):
    for fld in tqdm(all_flds):
        
        #ncols = len(im_)
        for icol in range(ncols-1):
            #icol=2
            tag = os.path.basename(fld)
            save_fl = save_folder+os.sep+fov.split('.')[0]+'--'+tag+'--col'+str(icol)+'__Xhfits.npz'
            if not os.path.exists(save_fl) or redo:
                im_ = read_im(fld+os.sep+fov)
                #print("Reading image")
                im__ = np.array(im_[icol],dtype=np.float32)
                
                im_n = norm_slice(im__,s=30)
                #print("Fitting image")
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
       r'\\192.168.0.3\bbfishdc10\DCBBL1_3_2_2023\MERFISH\H14_MER_set1',
       r'\\192.168.0.3\bbfishdc10_2\DCBBL1_3_2_2023\H15_MER_set1',
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
        print(len(all_flds),all_flds)
        compute_fits(save_folder,fov,all_flds,redo=False)
        print("Computing drift on: "+str(fov))
        compute_drift(save_folder,fov,all_flds,set_,redo=False)
    return set_ifov
if __name__ == '__main__':
    # start 4 worker processes
    items = [(set_,ifov)for set_ in ['_set1','_set2','_set3','_set4']
                        for ifov in range(1000)]
                        
    #main_f(['_set1',3])
    if True:
        with Pool(processes=35) as pool:
            print('starting pool')
            result = pool.map(main_f, items)

    