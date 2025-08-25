function vol2r = vol_resample_amd(vol2, vol1, M_vol1_to_vol2, varargin)
% Resample an image volume after a linear transformation. 
%
% function vol2r = vol_resample(vol2, vol1, M_vol1_to_vol2, [interpm], [padding], [bclamp])
%
%       vol2    : input image volume
%       vol1    : output image volume
%       M_vol1_to_vol2    : Reslice matrix from output volume to input volume
%       interpm : 0: : Nearest 1:Linear(default) 2:cubic 
%                 3: Key's spline 4: Cubic spline. 5: Hamming_Sinc
%       padding : half width of number of points used in the interpolation,
%                  only for 4 and 5. For example, 3 means it would use
%                  6*6*6 points for interpolation.
%       bclamp : set negative value to zero default = false
%   vol2r = vol_resample(vol2, vol1, M_vol1_to_vol2, [interpm]])
%
% NOTE 1: vol1's intensities are not used -- this volume simply
% defines the number of voxels and the "orientation" (e.g. coronal
% versus sagittal) of the output volume vol2r.
%
% NOTE 2: vol2r and vol1 will be in perfect register, assuming
% M_vol1_to_vol2 correctly describes how a point in vol1's LPH scanner
% coordinates maps to a point in vol2's LPH scanner coordinates.
%
% AMD: modified to handle n-d volumes

interpm = 1;
if nargin >= 4
  interpm = varargin{1};
end


padding=3;
if nargin >= 5
  padding = varargin{2};
end

bclamp=false;
if nargin >= 6
  bclamp = varargin{3};
end

if (interpm==0) % Nearest Neighbor padding=0
    padding=0;
elseif (interpm==1) % Lineae  padding=1
    padding=1;
elseif (interpm==2) % cubic and key's padding =2
    padding=2;
elseif (interpm==3)
    padding=2;
end

if iscell(vol2)
  vol2r = cell(size(vol2)); 
  for voli = 1:length(colvec(vol2))
    vol2r{voli} = vol1;
    if (isfield(vol2{voli}, 'dcminfo'))
      vol2r{voli}.dcminfo=vol2{voli}.dcminfo;
    else
      vol2r{voli}.dcminfo=[];
    end
    Mvxl2vxl=inv(vol2{voli}.Mvxl2lph)*M_vol1_to_vol2*vol2r{voli}.Mvxl2lph;
    vol2r{voli}.imgs=resliceMEX(double(vol2{voli}.imgs), size(vol2{voli}.imgs), Mvxl2vxl, size(vol2r{voli}.imgs), interpm, padding); % Added to conversion to double, to avoid crash -- should make resliceMEX handle single
    if (bclamp)
      ind=find(vol2r{voli}.imgs<0);
      vol2r{voli}.imgs(ind)=0;
    end
    [vol2r{voli}.maxI vol2r{voli}.minI] = maxmin(vol2r{voli}.imgs);
  end
else
  dims1 = size(vol1.imgs); 
  dims2 = size(vol2.imgs); if length(dims2)<4, dims2 = [dims2 1]; end
  nvols = prod(dims2(4:end));
  vol2r=vol1;
  if (isfield(vol2, 'dcminfo'))
      vol2r.dcminfo=vol2.dcminfo;
  else
      vol2r.dcminfo=[];
  end
  Mvxl2vxl=inv(vol2.Mvxl2lph)*M_vol1_to_vol2*vol2r.Mvxl2lph;
  vol2r_tmp = zeros([dims1(1:3) nvols]);
  for vi = 1:nvols
    vol2r_tmp(:,:,:,vi)=resliceMEX(double(vol2.imgs(:,:,:,vi)), dims2(1:3), Mvxl2vxl, dims1(1:3), interpm, padding); % Added to conversion to double, to avoid crash -- should make resliceMEX handle single
  end
  if length(dims2)>4
    vol2r_tmp = reshape(vol2r_tmp,[dims1(1:3) dims2(4:end)]);
  end
  vol2r.imgs=vol2r_tmp;
  if (bclamp)
      ind=find(vol2r.imgs<0);
      vol2r.imgs(ind)=0;
  end
  [vol2r.maxI vol2r.minI] = maxmin(vol2r.imgs);
end





