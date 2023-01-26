#Bogdan Bintu
#Copyright Presidents and Fellows of Harvard College, 2017.

import numpy as np
import glob,os
import matplotlib.pylab as plt
###These are aligment tools for either 3d/2d images or point sets

##3d images
def grab_block_im(im,center,block_sizes,return_coords=False):
    """Given an n-dim image <im>, a position <center> and a list of sizez <block_sizes>,
    retuns a block of the image of size <block_sizes> from <im> centered at <center>"""
    dims = im.shape
    coords = []
    def in_dim(c,dim):
        c_ = c
        if c_<0: c_=0
        if c_>dim: c_=dim
        return c_
    for c,block,dim in zip(center,block_sizes,dims):
        block_ = int(block/2)
        c=int(c)
        c_min,c_max = in_dim(c-block_,dim),in_dim(c+block-block_,dim)
        coords.append((c_min,c_max))
    slices = [slice(cm,cM) for cm,cM in coords]+[Ellipsis]
    if return_coords:
        return im[slices],coords
    return im[slices]
def grab_block(im,center,block_sizes,return_coords=False):
    """Given an n-dim image <im>, a position <center> and a list of sizez <block_sizes>,
    retuns a block of the image of size <block_sizes> from <im> centered at <center>"""
    if type(im) in [np.core.memmap,np.ndarray]:
        return grab_block_im(im,center,block_sizes,return_coords=return_coords)
        
    dims = im.shape
    coords = []
    def in_dim(c,dim):
        c_ = c
        if c_<0: c_=0
        if c_>dim: c_=dim
        return c_
    for c,block,dim in zip(center,block_sizes,dims):
        block_ = int(block/2)
        c=int(c)
        c_min,c_max = in_dim(c-block_,dim),in_dim(c+block-block_,dim)
        coords.append((c_min,c_max))
    ((zmin,zmax),(ymin,ymax),(xmin,xmax)) = coords
    im_sm = im.get_im(ind_col=None,minx=xmin,maxx=xmax,miny=ymin,maxy=ymax)[zmin:zmax].swapaxes(1,2)
    if return_coords:
        return im_sm,coords
    return im_sm
def fftalign(im1,im2,dm=100,plt_val=False):
    """
    Inputs: 2 3d images (numpy arrays) im1, im2 and a maximum displacement dm.
    This computes the cross-cor between im1 and im2 using numpy's fftconvolve (fast) and determines the maximum value within dm in all dimensions of 0.
    It returns the displacement vector <zxy_t> between im2 and im1.
    To be used with: translate(im2,zxy_t)
    if plt_val is True then a plot of the correlation matrix for the last 2 dimensions (generally x,y) is shown, together with the displacement 
    """    
    from scipy.signal import fftconvolve
    sh = np.array(im2.shape)
    dim1,dim2 = np.max([sh-dm,sh*0],0),sh+dm
    im2_=np.array(im2[dim1[0]:dim2[0],dim1[1]:dim2[1],dim1[2]:dim2[2]][::-1,::-1,::-1],dtype=float)
    im2_-=np.mean(im2_)
    im2_/=np.std(im2_)
    sh = np.array(im1.shape)
    dim1,dim2 = np.max([sh-dm,sh*0],0),sh+dm
    im1_=np.array(im1[dim1[0]:dim2[0],dim1[1]:dim2[1],dim1[2]:dim2[2]],dtype=float)
    im1_-=np.mean(im1_)
    im1_/=np.std(im1_)
    im_cor = fftconvolve(im1_,im2_, mode='full')
       
    xyz = np.unravel_index(np.argmax(im_cor), im_cor.shape)
    if np.sum(im_cor>0)>0:
        im_cor[im_cor==0]=np.min(im_cor[im_cor>0])
    else:
        im_cor[im_cor==0]=0
    if plt_val:
        plt.figure()
        x,y=xyz[-2:]
        im_cor_2d = np.array(im_cor)
        while len(im_cor_2d.shape)>2:
            im_cor_2d = np.max(im_cor_2d,0)
        plt.plot([y],[x],'ko')
        plt.imshow(im_cor_2d,interpolation='nearest')
        plt.show()
    xyz=np.round(-np.array(im_cor.shape)/2.+xyz).astype(int)
    return xyz

def translate(im,trans):
    """translates 3d images"""
    shape_ = im.shape
    zmax=shape_[0]
    xmax=shape_[1]
    ymax=shape_[2]
    zmin,xmin,ymin=0,0,0
    trans_=np.array(np.round(trans),dtype=int)
    zmin-=trans_[0]
    zmax-=trans_[0]
    xmin-=trans_[1]
    xmax-=trans_[1]
    ymin-=trans_[2]
    ymax-=trans_[2]
    im_base_0 = np.zeros([zmax-zmin,xmax-xmin,ymax-ymin])
    im_zmin = min(max(zmin,0),shape_[0])
    im_zmax = min(max(zmax,0),shape_[0])
    im_xmin = min(max(xmin,0),shape_[1])
    im_xmax = min(max(xmax,0),shape_[1])
    im_ymin = min(max(ymin,0),shape_[2])
    im_ymax = min(max(ymax,0),shape_[2])
    im_base_0[(im_zmin-zmin):(im_zmax-zmin),(im_xmin-xmin):(im_xmax-xmin),(im_ymin-ymin):(im_ymax-ymin)]=im[im_zmin:im_zmax,im_xmin:im_xmax,im_ymin:im_ymax]
    return im_base_0

def align_and_shift(im,imref,imtarget=None,dm=100,plt_val=False,verbose=False):
    """Performs val = fftalign(imref,im,dm,plt_val)
    and imtrans=translate(imtarget,-val)
    imtarget can be differnet than im. If None, it defaults to im
    Returns imtrans
    """

    val = fftalign(imref,im,dm=dm)
    if verbose:
        print(val)
    if imtarget is None:
        imtarget = im
    imtrans=translate(imtarget,-val)
    if plt_val:
        f, (ax1, ax2) = plt.subplots(1, 2, sharex=True,sharey=True)
        ax1.imshow(np.max(imref,0),interpolation='nearest',cmap='gray')
        ax1.set_title('Original')
        ax2.imshow(np.max(translate(im,-val),0),interpolation='nearest',cmap='gray')
        ax2.set_title('Translated')
        plt.show()
    return imtrans
    
##2d images
import cv2
def blurnorm2d(im,gb):
    """Normalize an input 2d image <im> by dividing by a cv2 gaussian filter of the image"""
    im_ = im.astype(np.float32)
    blurred = cv2.blur(im_,(gb,gb))
    return im_/blurred
def fft3d_from2d(im1,im2,gb=5,max_disp=150):
    """Given a refence 3d image <im1> and a target image <im2> 
    this max-projects along the first (z) axis and finds the best tx,ty using fftalign_2d.
    Then it trims and max-projects along the last (y) axis and finds tz.
    Before applying fftalignment we normalize the images using blurnorm2d for stability."""
    im1_ = blurnorm2d(np.max(im1,0),gb)
    im2_ = blurnorm2d(np.max(im2,0),gb)
    tx,ty = fftalign_2d(im1_, im2_, center=[0, 0], max_disp=max_disp, plt_val=False)
    sx,sy = im1_.shape
    im1_t = blurnorm2d(np.max(im1[:,max(tx,0):sx+tx,max(ty,0):sy+ty],axis=-1),gb)
    im2_t = blurnorm2d(np.max(im2[:,max(-tx,0):sx-tx,max(-ty,0):sy-ty],axis=-1),gb)
    tz,_ = fftalign_2d(im1_t, im2_t, center=[0, 0], max_disp=max_disp, plt_val=False)
    return np.array([tz,tx,ty])
def fftalign_2d(im1,im2,center=[0,0],max_disp=50,return_cor_max=False,plt_val=False):
    """
    Inputs: 2 2D images <im1>, <im2>, the expected displacement <center>, the maximum displacement <max_disp> around the expected vector.
    This computes the cross-cor between im1 and im2 using fftconvolve (fast) and determines the maximum
    """    
    from scipy.signal import fftconvolve
    im2_=np.array(im2[::-1,::-1],dtype=float)
    im2_-=np.mean(im2_)
    im2_/=np.std(im2_)
    im1_=np.array(im1,dtype=float)
    im1_-=np.mean(im1_)
    im1_/=np.std(im1_)
    im_cor = fftconvolve(im1_,im2_, mode='full')

    sx_cor,sy_cor = im_cor.shape
    center_ = np.array(center)+np.array([sx_cor,sy_cor])/2.
    
    x_min = int(min(max(center_[0]-max_disp,0),sx_cor))
    x_max = int(min(max(center_[0]+max_disp,0),sx_cor))
    y_min = int(min(max(center_[1]-max_disp,0),sy_cor))
    y_max = int(min(max(center_[1]+max_disp,0),sy_cor))
    
    im_cor0=np.zeros_like(im_cor)
    im_cor0[x_min:x_max,y_min:y_max]=1
    im_cor = im_cor*im_cor0
       
    y, x = np.unravel_index(np.argmax(im_cor), im_cor.shape)
    if np.sum(im_cor>0)>0:
        im_cor[im_cor==0]=np.min(im_cor[im_cor>0])
    else:
        im_cor[im_cor==0]=0
    if plt_val:
        plt.figure()
        plt.plot([x],[y],'k+')
        plt.imshow(im_cor,interpolation='nearest')
        plt.xlim(im_cor.shape[0]/2-max_disp,im_cor.shape[0]/2+max_disp)
        plt.ylim(im_cor.shape[1]/2-max_disp,im_cor.shape[1]/2+max_disp)
        plt.show()
    xt,yt=(-(np.array(im_cor.shape)-1)/2.+[y,x]).astype(int)
    if return_cor_max:
        return (xt,yt),np.max(im_cor)
    return xt,yt
##point clouds

def translation_aling_pts(cents_fix,cents_target,cutoff=2.,xyz_res=1,plt_val=False,return_pts=False):
    """
    This checks all pairs of points in cents_target for counterparts of same distance (+/- cutoff) in cents_fix
    and adds them as posibilities. Then uses multi-dimensional histogram across txyz with resolution xyz_res.
    Then it finds nearest neighbours and returns the median txyz_b within resolution.
    
    
    Returns the translation vector such that cents_target + txyz_f ->cents_fix, 
    if return_pts it also returns 
    
    Warning: This is intrinsically O(N^4) which can get quite big. Expect ~1s runtime for ~100 points.
    Although slow, this procere is robust.
    """
    #Lazy imports
    
    from itertools import combinations
    from scipy.spatial.distance import pdist
    from scipy.spatial.distance import cdist
    
    
    if len(cents_fix)==0 or len(cents_target)==0:
        if return_pts:
            return np.array([np.inf]*3),[],[]
        return np.array([np.inf]*3)

    cents = np.array(cents_fix)
    cents_target = np.array(cents_target)
    dists_target = pdist(cents_target)
    dists = pdist(cents_fix)
    all_pairs = np.array(list(combinations(range(len(cents)),2)))
    all_pairs_target = np.array(list(combinations(range(len(cents_target)),2)))
    #inds_all = np.arange(len(dists))
    txyzs=[]
    for ind_target in range(len(dists_target)):
        keep_cands = np.abs(dists-dists_target[ind_target])<cutoff
        good_pairs = all_pairs[keep_cands][:]
        if len(good_pairs)>0:
            p1 = cents[good_pairs[:,0]]
            p2 = cents[good_pairs[:,1]]
            p1T = cents_target[all_pairs_target[ind_target,0]]
            p2T = cents_target[all_pairs_target[ind_target,1]]
            txyzs.extend(p1[:]-[p1T])
            txyzs.extend(p1[:]-[p2T])
    if len(txyzs)==0:
        if return_pts:
            return np.array([np.inf]*3),[],[]
        return np.array([np.inf]*3)
    bin_txyz = np.array((np.max(txyzs,axis=0)-np.min(txyzs,axis=0))/float(xyz_res),dtype=int)+1
    #print bin_txyz
    hst_res = np.histogramdd(np.array(txyzs),bins=bin_txyz)
    #print "Finished Histogram"
    if plt_val:
        plt.figure()
        plt.imshow(np.max(hst_res[0],0),interpolation='nearest')
        plt.show()
    #print "Finding best guess value"
    ibest = np.unravel_index(np.argmax(hst_res[0]),hst_res[0].shape)
    txyz_f = [hst[ib]for hst,ib in zip(hst_res[1],ibest)]
    txyz_f = np.array(txyz_f)
    #print "Pick closest matches"
    inds_closestT = np.argmin(cdist(cents,cents_target + txyz_f),axis=1)
    inds_closestF=np.arange(len(inds_closestT))
    #print "Filter closest matches"
    keep = np.sqrt(np.sum((cents_target[inds_closestT]+ txyz_f-cents[inds_closestF])**2,axis=-1))<2*xyz_res
    inds_closestT=inds_closestT[keep]
    inds_closestF=inds_closestF[keep]
    #print "Compute and return median"
    txyz_b = np.median(cents_target[inds_closestT]-cents[inds_closestF],axis=0)
    if plt_val:
        plt.figure()
        plt.plot(cents[inds_closestF].T[0],cents[inds_closestF].T[1],'go')
        plt.plot(cents_target[inds_closestT].T[0]-txyz_b[0],cents_target[inds_closestT].T[1]-txyz_b[1],'ro')
        plt.figure()
        dists = np.sqrt(np.sum((cents_target[inds_closestT]-cents[inds_closestF])**2,axis=-1))
        plt.hist(dists)
        plt.show()
    if return_pts:
        return txyz_b,cents[inds_closestF],cents_target[inds_closestT]
    return txyz_b