function vol_ctx = niftiread_ctx(fname_nii)

[vol M] = niftiread_amd(fname_nii);
vol_ctx = ctx_mgh2ctx(vol,M);

